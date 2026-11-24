#!/usr/bin/env python3
"""A12 passo 2: termo de vizinhos ausentes nos 29 orbitais (s3p3d2f1), só do volume (TeBulk_fo2).
ΔH_μν(i) = −Σ_aus ⟨μ|f_VNA(|r − R_k|)|ν⟩ + ⟨μ|δV_xc|ν⟩, δV_xc = V_xc^PBE[ρ_vol − Σ_aus ρ_at] − V_xc^PBE[ρ_vol].
f_VNA: ajuste em espaço real do A11 r2 (a11_r2_fvna.npz, R_c 6 bohr); ρ_vol = tden.cube (fo2); ρ_at como no A8/A11.
Malha: a do cubo do OpenMX (igual ao A11), numa caixa local em torno dos centros (valores do volume por índice periódico),
margem de 3 pontos para as diferenças centrais (grad + div); nada de imagens periódicas para os ausentes.
Orbitais: ângulos de pao_eval (mod_pao_basis.f90, ordem/sinais do OpenMX: l, ζ, m); radial por spline cúbica do .pao
(conferido contra o próprio pao_eval compilado, paochk/paoeval_drv).
Uso: a12_termo.py check  → confere com o A11 (blocos s, p, átomos 0/1/30/32 do N11)
     a12_termo.py film   → termo para os 12 átomos centrais do filme L8 (colunas 4..7) + hoppings da face; grava .sitepot"""
import sys, subprocess, numpy as np
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
# raízes configuráveis (padrão: layout original da campanha)
LFI_TESE = _os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/")
PAO = _os.environ.get("LFI_PAO", _os.path.join(_HERE, "..", "..", "test", "data", "pao", "Te7.0.pao"))
D = LFI_TESE
T = D + "campo_semiinf/termo/"
sys.path.insert(0, _HERE)
from a11_pbe import GridFD, vxc
import a8lib
HA = 27.211386245988; B2A = 0.529177210903
RORB = 7.0

def readcube(f):
    L = open(f).read().split("\n")
    nat = int(L[2].split()[0]); org = np.array(L[2].split()[1:4], float)
    n = np.array([int(L[3 + k].split()[0]) for k in range(3)])
    ax = np.array([L[3 + k].split()[1:4] for k in range(3)], float)
    tok = " ".join(L[6 + nat:]).split()
    return np.array(tok, float).reshape(n), org, n, ax
rho, org, ncub, ax = readcube(_os.environ.get("LFI_TDEN", "/home/nathanael/calc/tese-Te/green_surface/TeBulk_fo2/TeBulk.tden.cube"))
dVol = abs(np.linalg.det(ax))
iaxT = np.linalg.inv(ax.T)
fz = np.load(D + "vizinhos/a11_r2_fvna.npz"); rr, fV = fz["r"], fz["f"]
xv, den = a8lib.atomic_density(pao=PAO, zval=16.0); frho = a8lib.make_rho_atom(xv, den)

# ---------- PAO (s3p3d2f1): o próprio pao_eval do LowFoldIn (paochk/paoeval_bin) ----------
import os, tempfile
def orbitals(rel):
    with tempfile.TemporaryDirectory() as td:
        fi, fo = os.path.join(td, "in.bin"), os.path.join(td, "out.bin")
        with open(fi, "wb") as f:
            np.array([len(rel)], np.int32).tofile(f); np.ascontiguousarray(rel, np.float64).tofile(f)
        subprocess.run([_os.environ.get("LFI_PAOEVAL", T + "paochk/paoeval_bin"), fi, fo, PAO], check=True)
        return np.fromfile(fo, np.float64).reshape(len(rel), 29).T
def check_pao():   # pz e s contra a fórmula do A11 (radial linear do .pao): ordem e sinais
    P = open(PAO).read().split("\n")
    def block(tag):
        i0 = [k for k, l in enumerate(P) if l.strip().startswith("<" + tag)][0]
        i1 = [k for k, l in enumerate(P) if l.strip().startswith(tag + ">")][0]
        return np.array([l.split() for l in P[i0 + 1:i1]], float)
    R0 = block("pseudo.atomic.orbitals.L=0"); R1 = block("pseudo.atomic.orbitals.L=1")
    d = np.random.default_rng(1).uniform(-5, 5, (400, 3)); r = np.linalg.norm(d, axis=1)
    ref = [np.interp(r, R0[:, 1], R0[:, 2]) / np.sqrt(4 * np.pi)] + [np.sqrt(3 / (4 * np.pi)) * np.interp(r, R1[:, 1], R1[:, 2]) * d[:, k] / r for k in range(3)]
    F = orbitals(d)
    return np.abs(F[[0, 3, 4, 5]] - np.array(ref)).max(), np.abs(F).max()

# ---------- caixa local ----------
def box(centers, Rb):
    """pontos da malha do cubo (bohr) a < Rb de algum centro, com caixa + margem 3; devolve rS, rhoS, máscara."""
    g = np.array([iaxT @ (c - org) for c in centers])
    ext = Rb * np.linalg.norm(np.linalg.inv(ax.T), axis=1)   # |Δi_k| ≤ Rb ‖linha k de ax^{-T}‖
    lo = np.floor(g.min(0) - ext).astype(int) - 3; hi = np.ceil(g.max(0) + ext).astype(int) + 3
    nb = hi - lo + 1
    I = np.stack(np.meshgrid(*[np.arange(lo[k], hi[k] + 1) for k in range(3)], indexing="ij"), -1).reshape(-1, 3)
    rS = org + I @ ax
    rhoS = rho[I[:, 0] % ncub[0], I[:, 1] % ncub[1], I[:, 2] % ncub[2]]
    inside = np.zeros(len(rS), bool)
    for c in centers: inside |= np.sum((rS - c) ** 2, 1) < Rb ** 2
    return rS, rhoS, inside, nb

def termo(Ri, mis, Rj=None):
    """blocos (eV) ⟨φ_i|V|φ_j⟩ (j = i se Rj é None): S, VNA, δV_xc PBE. Ri, Rj, mis em bohr."""
    cen = [Ri] if Rj is None else [Ri, Rj]
    rS, rhoS, inside, nb = box(cen, RORB)
    rsub = rhoS.copy(); pot = np.zeros(len(rS))
    for k in mis:
        dk = np.linalg.norm(rS - k, axis=1)
        m = dk < a8lib.RCUT_ATOM; rsub[m] -= frho(dk[m])
        m = dk < rr[-1]; pot[m] -= np.interp(dk[m], rr, fV)
    g = GridFD(ax * nb[:, None], nb)
    dvx = (vxc(rsub.reshape(nb), g, True) - vxc(rhoS.reshape(nb), g, True)).ravel()
    sel = inside
    phi = orbitals(rS[sel] - Ri); phj = phi if Rj is None else orbitals(rS[sel] - Rj)
    S = phi @ phj.T * dVol
    Vna = (phi * pot[sel]) @ phj.T * dVol * HA
    Vxc = (phi * dvx[sel]) @ phj.T * dVol * HA
    return S, Vna, Vxc, int(np.sum(rsub[sel] < 1e-6))

s1p1 = [0, 3, 4, 5]
def metrics(dH, S):
    cc = np.sum(S[np.ix_(s1p1, s1p1)] * dH[np.ix_(s1p1, s1p1)]) / np.sum(S[np.ix_(s1p1, s1p1)] ** 2)
    pb = dH[3:6, 3:6]; ap = pb - np.trace(pb) / 3 * np.eye(3); w, U = np.linalg.eigh(ap)
    k = np.argmax(np.abs(w)); oth = [j for j in range(3) if j != k]
    return cc, w[k] - w[oth].mean(), U[:, k], dH[0, 3:6], np.linalg.norm(dH[np.ix_([0, 1, 2], list(range(3, 12)))])
def fmt(M, S, ref=None):
    cc, sp_, axv, vec, nsp = metrics(M, S)
    s = f"{cc:+.4f} {sp_:+.4f} | {np.linalg.norm(vec):.4f} | {nsp:.4f}"
    if ref is not None:
        md = metrics(ref, S); s += f" | cos(eixo) {abs(axv @ md[2]):.2f} cos(v) {vec @ md[3] / (np.linalg.norm(vec) * np.linalg.norm(md[3]) + 1e-30):+.2f}"
    return s

if sys.argv[1] == "check":
    e, fm = check_pao(); print(f"# pao: max|φ_py − φ_pao_eval| = {e:.2e} (max|φ| {fm:.2f}) em 400 pontos")
    dat = np.load(D + "vizinhos/a11_blocos.npz"); SP = list(range(12))
    print("# át n | fonte | escalar σ−π | |v| | ‖sp‖ (eV)   — comparar com a11_r2_saida.txt")
    for i in (0, 1, 30, 32):
        dk = dat[f"mis{i}"] / B2A; Ri = (dat["pos"][i] - dat["t"]) / B2A
        S, Vna, Vxc, neg = termo(Ri, [Ri + d for d in dk])
        Sd = dat[f"S{i}"][np.ix_(SP, SP)]; dHd = (dat[f"Hs{i}"] - dat[f"Hb{i}"])[np.ix_(SP, SP)]
        print(f"#   át {i}: max|S_malha − S_scfout| (29×29) {np.abs(S[:12, :12] - Sd).max():.1e}; ρ_sub < 1e-6 em {neg} pts")
        for k, M in (("VNA", Vna), ("δV_xc PBE", Vxc), ("TOTAL", Vna + Vxc)):
            print(f"{i:3d} {len(dk)} | {k:10s} | {fmt(M[:12, :12], Sd, dHd)}")
    sys.exit()

# ---------- filme L8 ----------
def readstruct(f):
    L = open(f).read().split("\n")
    lat = np.array([[float(x) for x in l.split("=")[1].split()] for l in L if "LatticeVectors" in l])
    i = [k for k, l in enumerate(L) if l.strip() == "&ATOMS"][0]; nat = int(L[i + 1])
    return lat, np.array([[float(x) for x in L[i + 2 + a].split()[1:4]] for a in range(nat)])
latb, posb = readstruct(T + "hs/TeBulk.structure.nml")
latf, posf = readstruct(T + "hs/TeFilm_x1_L8.structure.nml")
# filme = sítios do volume com n2 = 0..7 (confere); ausentes = n2 ≥ 8
for a in range(24):
    assert np.abs(posf[a] - (posb[a % 3] + (a // 3) * latb[1])).max() < 1e-6
RMAX = RORB + max(6.0, a8lib.RCUT_ATOM) + 0.5   # 14,5 bohr: orbital 7 + max(R_c VNA 6, ρ_at 7)
cent = list(range(12, 24))
def missing(Ri_A):
    out = []
    for n1 in range(-5, 6):
        for n2 in range(8, 14):
            for n3 in range(-4, 5):
                for b in range(3):
                    R = posb[b] + n1 * latb[0] + n2 * latb[1] + n3 * latb[2]
                    if np.linalg.norm(R - Ri_A) / B2A < RMAX: out.append(R / B2A)
    return out
blocks = {}
print(f"# filme L8, átomos centrais 12..23 (colunas 4..7); ausentes: sítios n2 ≥ 8 a < {RMAX:.1f} bohr")
print("# át col nmis d_min(Å) | escalar σ−π | |v| | ‖sp‖ (s,p; eV) | máx|ΔH| 29×29 | máx|ΔH| d,f | máx|S−I|")
for a in cent:
    Ri = posf[a] / B2A; mis = missing(posf[a])
    if not mis: continue
    S, Vna, Vxc, neg = termo(Ri, mis)
    dH = Vna + Vxc; dH = 0.5 * (dH + dH.T); blocks[a] = dH
    dmin = min(np.linalg.norm(m - Ri) for m in mis) * B2A
    df = np.abs(dH[12:, :]).max()
    print(f"{a:3d} {a // 3} {len(mis):3d} {dmin:5.2f} | {fmt(dH[:12, :12], S[:12, :12])} | {np.abs(dH).max():.4f} | {df:.4f} | {np.abs(S - np.eye(29)).max():.1e}")
np.savez(T + "a12_blocos.npz", **{f"dH{a}": blocks[a] for a in blocks}, pos=posf)
# hoppings entre vizinhos da face (coluna 7 e 7–6), d < 3,2 Å
print("# fora da diagonal: pares vizinhos (d < 3,2 Å, mesma célula em x/z ou imagem), máx|⟨μ_i|V|ν_j⟩| × máx on-site")
mx = 0.0
for a in range(18, 24):
    for b in range(a + 1, 24):
        for n1 in (-1, 0, 1):
            for n3 in (-1, 0, 1):
                Rb = posf[b] + n1 * latf[0] + n3 * latf[2]
                if np.linalg.norm(Rb - posf[a]) > 3.2: continue
                Ri, Rj = posf[a] / B2A, Rb / B2A
                mis = [m for m in missing(posf[a])] + [m for m in missing(Rb) if all(np.linalg.norm(m - q) > 1e-6 for q in missing(posf[a]))]
                S, Vna, Vxc, _ = termo(Ri, mis, Rj)
                h = np.abs(Vna + Vxc).max(); mx = max(mx, h)
                print(f"#   {a}-{b} ({n1},{n3}) d = {np.linalg.norm(Rb - posf[a]):.2f} Å: máx|ΔH_ij| = {h:.4f} eV; on-site máx {max(np.abs(blocks.get(a, np.zeros(1))).max(), np.abs(blocks.get(b, np.zeros(1))).max()):.4f}")
# arquivo LFSITEPOT (bloco por átomo, eV)
with open(T + "termo_aus.sitepot", "w") as f:
    f.write("LFSITEPOT 1\n# A12: -sum_aus <VNA_at> + <dV_xc PBE> (TeBulk_fo2), s3p3d2f1, eV; posições do filme L8 (Å)\n")
    f.write(f"norb 29\nnentry {len(blocks)}\n")
    for a, M in blocks.items():
        f.write(f"entry 0 {a + 1} {posf[a][0]:.10f} {posf[a][1]:.10f} {posf[a][2]:.10f} block\n")
        for row in M: f.write(" ".join(f"{v:.10f}" for v in row) + "\n")
with open(T + "zero.sitepot", "w") as f:
    f.write("LFSITEPOT 1\n# bloco zero (teste bit a bit)\n")
    f.write(f"norb 29\nnentry {len(blocks)}\n")
    for a in blocks:
        f.write(f"entry 0 {a + 1} {posf[a][0]:.10f} {posf[a][1]:.10f} {posf[a][2]:.10f} block\n")
        for row in range(29): f.write(" ".join("0.0" for _ in range(29)) + "\n")
print(f"# máx|ΔH_ij| (pares) = {mx:.4f} eV")
print("# gravados termo_aus.sitepot, zero.sitepot, a12_blocos.npz")
