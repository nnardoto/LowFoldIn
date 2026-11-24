#!/usr/bin/env python3
"""A18: barreira suave tipo Fermi–Dirac no vácuo, V_b(y) = V0 / (1 + exp(−(y − y0)/w)), y0 = y_topo + d (y normal, vácuo em +y;
y_topo = átomo central mais externo do filme L8). Blocos on-site ΔH_μν(i) = ⟨μ_i|V_b|ν_i⟩ (29 orbitais, s3p3d2f1) com o avaliador de
PAO do A12 (paochk/paoeval_bin = pao_eval do LowFoldIn), malha cartesiana de passo h (bohr) na esfera de 7 bohr; somados ao termo
VNA+xc do A12 (termo_aus.sitepot) → LFSITEPOT 'block'. Só entram átomos com máx|bloco| > 1 meV. Pares i≠j vizinhos (d < 3,2 Å):
máx|⟨μ_i|V_b|ν_j⟩| medido e comparado ao on-site (só reportado; o formato aceita só on-site).
uso: a18_barreira.py V0_eV d_A w_A out.sitepot [h_bohr] [pares]   |   a18_barreira.py conv  (convergência em h)"""
import sys, os, subprocess, tempfile, numpy as np
B2A = 0.529177210903
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
# raízes configuráveis (padrão: layout original da campanha)
LFI_TESE = _os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/")
PAO = _os.environ.get("LFI_PAO", _os.path.join(_HERE, "..", "..", "test", "data", "pao", "Te7.0.pao"))
T = LFI_TESE + "campo_semiinf/termo/"
PB = _os.environ.get("LFI_PAOEVAL", T + "paochk/paoeval_bin"); RC = 7.0
def orbitals(rel):
    with tempfile.TemporaryDirectory() as td:
        fi, fo = os.path.join(td, "in.bin"), os.path.join(td, "out.bin")
        with open(fi, "wb") as f:
            np.array([len(rel)], np.int32).tofile(f); np.ascontiguousarray(rel, np.float64).tofile(f)
        subprocess.run([PB, fi, fo, PAO], check=True)
        return np.fromfile(fo, np.float64).reshape(len(rel), 29).T
def readstruct(f):
    L = open(f).read().split("\n")
    lat = np.array([[float(x) for x in l.split("=")[1].split()] for l in L if "LatticeVectors" in l])
    i = [k for k, l in enumerate(L) if l.strip() == "&ATOMS"][0]; nat = int(L[i + 1])
    return lat, np.array([[float(x) for x in L[i + 2 + a].split()[1:4]] for a in range(nat)])
latf, posf = readstruct(T + "hs/TeFilm_x1_L8.structure.nml")
cent = list(range(12, 24)); ytop = posf[cent, 1].max()
_G = {}
def grid(h):
    if h not in _G:
        x = np.arange(-RC, RC + h / 2, h); X, Y, Z = np.meshgrid(x, x, x, indexing="ij"); s = X**2 + Y**2 + Z**2 <= RC**2
        P = np.stack([X[s], Y[s], Z[s]], 1); _G[h] = (P, orbitals(P), h**3)
    return _G[h]
def vb(yA, V0, d, w):
    return V0 / (1 + np.exp(-np.clip((yA - (ytop + d)) / w, -700, 700)))
def onsite(i, V0, d, w, h):
    P, phi, h3 = grid(h); v = vb(posf[i, 1] + P[:, 1] * B2A, V0, d, w)
    M = (phi * v) @ phi.T * h3; return 0.5 * (M + M.T), phi @ phi.T * h3
def read_sitepot(f):
    L = open(f).read().split("\n"); out = {}; k = 0
    while k < len(L):
        if L[k].startswith("entry"):
            a = int(L[k].split()[2]) - 1
            out[a] = np.array([L[k + 1 + r].split() for r in range(29)], float); k += 30
        else: k += 1
    return out
if sys.argv[1] == "conv":
    for h in (0.20, 0.15, 0.12):
        Ms = [onsite(i, 1.0, 0.5, 0.5, h)[0] for i in cent]
        S = onsite(cent[0], 0, 0, 0.5, h)[1]
        if h == 0.20: M0 = Ms
        print(f"h = {h:.2f} bohr: máx|S − I| (diag) {np.abs(np.diag(S) - 1).max():.1e}; máx|ΔM vs h=0,20| {max(np.abs(a - b).max() for a, b in zip(Ms, M0)):.2e} eV;"
              f" máx|M| {max(np.abs(a).max() for a in Ms):.4f} eV (V0 = 1, d = 0,5, w = 0,5)")
    sys.exit()
V0, d, w, out = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
h = float(sys.argv[5]) if len(sys.argv) > 5 else 0.15
aus = read_sitepot(T + "termo_aus.sitepot")
print(f"# A18: V0 = {V0} eV, d = {d} Å, w = {w} Å, y0 = y_topo + d = {ytop:.4f} + {d} Å; malha h = {h} bohr")
print("# át y(Å) y−y_topo | máx|ΔH_b| (eV) | ⟨s1|V|s1⟩ ⟨p1y|V|p1y⟩ ⟨s1|V|p1y⟩ | entra")
blocks = {}
for i in cent:
    M, S = onsite(i, V0, d, w, h); mx = np.abs(M).max(); ok = mx > 1e-3
    print(f"{i + 1:3d} {posf[i, 1]:8.3f} {posf[i, 1] - ytop:+7.3f} | {mx:.4f} | {M[0, 0]:+.4f} {M[4, 4]:+.4f} {M[0, 4]:+.4f} | {'sim' if ok else 'não'}")
    if ok: blocks[i] = M
if len(sys.argv) > 6:
    P, phi, h3 = grid(h)
    print("# pares i≠j vizinhos (d < 3,2 Å) com algum on-site > 1 meV: máx|⟨μ_i|V_b|ν_j⟩| e razão ao máx on-site de i, j")
    rmax = 0
    for a in cent:
        for b in cent:
            if b <= a: continue
            for n1 in (-1, 0, 1):
                for n3 in (-1, 0, 1):
                    R = posf[b] + n1 * latf[0] + n3 * latf[2] - posf[a]; dd = np.linalg.norm(R)
                    if dd > 3.2 or (a not in blocks and b not in blocks): continue
                    phj = orbitals(P - R / B2A); v = vb(posf[a, 1] + P[:, 1] * B2A, V0, d, w)
                    Vij = (phi * v) @ phj.T * h3; on = max(np.abs(blocks[k]).max() for k in (a, b) if k in blocks)
                    r = np.abs(Vij).max() / on; rmax = max(rmax, r)
                    print(f"#   {a + 1}-{b + 1} ({n1},{n3}) d = {dd:.2f} Å: máx|V_ij| = {np.abs(Vij).max():.4f} eV = {100 * r:.0f} % do on-site máx {on:.4f}")
    print(f"# razão máxima par/on-site = {100 * rmax:.0f} %")
tot = {a: aus.get(a, 0) + blocks.get(a, 0) for a in sorted(set(aus) | set(blocks))}
with open(out, "w") as f:
    f.write("LFSITEPOT 1\n# A18: termo VNA+xc do A12 (termo_aus.sitepot) + barreira FD V0/(1+exp(-(y-y0)/w)), " +
            f"V0 = {V0} eV, d = {d} A, w = {w} A, y0 = {ytop + d:.6f} A; s3p3d2f1, eV; filme L8 (A)\n")
    f.write(f"norb 29\nnentry {len(tot)}\n")
    for a, M in tot.items():
        f.write(f"entry 0 {a + 1} {posf[a][0]:.10f} {posf[a][1]:.10f} {posf[a][2]:.10f} block\n")
        for row in M: f.write(" ".join(f"{x:.10f}" for x in row) + "\n")
print(f"escrito {out}: {len(tot)} entradas ({len(blocks)} com barreira, {len(aus)} do A12)")
