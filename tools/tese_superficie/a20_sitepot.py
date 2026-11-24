#!/usr/bin/env python3
"""A20: LFSITEPOT 'block' = termo VNA+xc do A12 + barreira FD do A18 (a18_barreira.py, V0, d, w) + L=2 uniaxial do par isolado
(A16/A18): c (ŷŷ − I/3) no bloco p1 (orbitais 4–6, px py pz) dos átomos 22 e 24 (face), ŷ = normal ao plano das duas ligações
covalentes (vizinhos < 3,2 Å, com imagens no plano da face). c = 0 → só barreira (controle). V0 = 0 → só VNA+xc (+ L=2).
Tensor opcional: '--tensor arq.npy' (dict átomo → 3×3 em eV) substitui c(ŷŷ − I/3) (rota só-volume da parte 2).
uso: a20_sitepot.py V0_eV d_A w_A c_eV out.sitepot [--tensor arq.npy]"""
import sys, os, subprocess, tempfile, numpy as np
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
# raízes configuráveis (padrão: layout original da campanha)
LFI_TESE = _os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/")
PAO = _os.environ.get("LFI_PAO", _os.path.join(_HERE, "..", "..", "test", "data", "pao", "Te7.0.pao"))
BAR = _os.path.join(_HERE, "a18_barreira.py")
T = LFI_TESE + "campo_semiinf/termo/"
def readstruct(f):
    L = open(f).read().split("\n")
    lat = np.array([[float(x) for x in l.split("=")[1].split()] for l in L if "LatticeVectors" in l])
    i = [k for k, l in enumerate(L) if l.strip() == "&ATOMS"][0]; nat = int(L[i + 1])
    return lat, np.array([[float(x) for x in L[i + 2 + a].split()[1:4]] for a in range(nat)])
latf, posf = readstruct(T + "hs/TeFilm_x1_L8.structure.nml")
def bond_normal(i, rcov=3.2):
    b = []
    for j in range(len(posf)):
        for n1 in (-1, 0, 1):
            for n3 in (-1, 0, 1):
                R = posf[j] + n1 * latf[0] + n3 * latf[2] - posf[i]; dd = np.linalg.norm(R)
                if 0.1 < dd < rcov: b.append(R)
    assert len(b) == 2, (i + 1, len(b))
    n = np.cross(b[0], b[1]); return n / np.linalg.norm(n), b
def l2_blocks(c):
    out = {}
    for a in (21, 23):
        y, _ = bond_normal(a); out[a] = c * (np.outer(y, y) - np.eye(3) / 3)
    return out
def read_sitepot(f):
    L = open(f).read().split("\n"); hdr = [l for l in L[:3] if l.startswith("#")]; out = {}; k = 0
    while k < len(L):
        if L[k].startswith("entry"):
            a = int(L[k].split()[2]) - 1
            out[a] = (L[k], np.array([L[k + 1 + r].split() for r in range(29)], float)); k += 30
        else: k += 1
    return hdr, out
if __name__ == "__main__":
    V0, d, w, c, out = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    ten = None
    if "--tensor" in sys.argv: ten = np.load(sys.argv[sys.argv.index("--tensor") + 1], allow_pickle=True).item()
    with tempfile.TemporaryDirectory() as td:
        tmp = os.path.join(td, "b.sitepot")
        if V0 != 0:
            subprocess.run([sys.executable, BAR, str(V0), str(d), str(w), tmp], check=True, stdout=subprocess.DEVNULL)
            hdr, E = read_sitepot(tmp)
        else:
            hdr, E = read_sitepot(T + "termo_aus.sitepot")
    add = ten if ten is not None else (l2_blocks(c) if c != 0 else {})
    for a, M3 in add.items():
        if a not in E:
            p = posf[a]; E[a] = (f"entry 0 {a + 1} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f} block", np.zeros((29, 29)))
        line, M = E[a]; M = M.copy(); M[3:6, 3:6] += M3; E[a] = (line, M)
    with open(out, "w") as f:
        f.write("LFSITEPOT 1\n")
        if ten is not None: f.write(f"# + L=2 tensor só-volume (A20, {sys.argv[sys.argv.index('--tensor') + 1]}) no bloco p1, átomos " + ", ".join(str(a + 1) for a in add) + "\n")
        elif c != 0: f.write(f"# + L=2 uniaxial par isolado (A16): {c:+.2f} eV (ŷŷ−I/3) no bloco p1, átomos 22 e 24\n")
        for h in hdr: f.write(h + "\n")
        f.write(f"norb 29\nnentry {len(E)}\n")
        for a in sorted(E):
            line, M = E[a]; f.write(line + "\n")
            for row in M: f.write(" ".join(f"{x:.10f}" for x in row) + "\n")
    print(f"escrito {out}: {len(E)} entradas; L=2 em {[a + 1 for a in add]}")
