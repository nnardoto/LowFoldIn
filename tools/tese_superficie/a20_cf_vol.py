#!/usr/bin/env python3
"""A20 parte 2, alternativa (i): campo cristalino L=2 do bloco p1 on-site do VOLUME (DFT TeBulk, Hb do a11_blocos.npz)
decomposto em contribuições uniaxiais de dois centros por vizinho, A_tl = Σ_n β_shell(n) (d̂_n d̂_n − I/3), camadas:
intra-cadeia (2 a 2,83 Å) e entre cadeias (4 a 3,75 Å) [+ 3ª camada opcional]; ajuste por mínimos quadrados nos 3 átomos do
volume (5 componentes sem traço cada). Previsão na face: ΔH_cf = −Σ_aus β (d̂d̂ − I/3); o que falta ALÉM do termo VNA+xc
(termo_aus.sitepot, que já remove a parte VNA+xc dos ausentes) é R = ΔH_cf − tl(termo_aus p1). c_pred = (3/2)⟨ŷ|R|ŷ⟩.
Conferência (só diagnóstico, não entra no ajuste): DFT N11 (Hs − Hb) sem traço nos átomos 30 (1 aus) e 32 (3 aus).
uso: a20_cf_vol.py"""
import sys, numpy as np
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
# raízes configuráveis (padrão: layout original da campanha)
LFI_TESE = _os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/")
PAO = _os.environ.get("LFI_PAO", _os.path.join(_HERE, "..", "..", "test", "data", "pao", "Te7.0.pao"))
sys.path.insert(0, _HERE)
from a20_sitepot import bond_normal, posf, latf, read_sitepot
D = LFI_TESE
d = np.load(D + "vizinhos/a11_blocos.npz"); btv, bpos = d["btv"], d["bpos"]
def tl(M): return M - np.trace(M) / 3 * np.eye(3)
def U(v): v = v / np.linalg.norm(v); return np.outer(v, v) - np.eye(3) / 3
def neigh(b, rmax=4.6):
    out = []
    for n1 in range(-2, 3):
        for n2 in range(-2, 3):
            for n3 in range(-2, 3):
                for j in range(3):
                    R = bpos[j] + n1 * btv[0] + n2 * btv[1] + n3 * btv[2] - bpos[b]; r = np.linalg.norm(R)
                    if 0.1 < r < rmax: out.append(R)
    return out
SH = [(0, 3.0), (3.0, 4.0), (4.0, 4.6)]
def basis(b, nsh):
    nb = neigh(b); return [sum((U(R) for R in nb if SH[s][0] < np.linalg.norm(R) < SH[s][1]), np.zeros((3, 3))) for s in range(nsh)]
# sublattice de cada átomo da fatia N11 (para achar o Hb correspondente)
def sub(i):
    r = d["pos"][i] - d["t"]
    for b in range(3):
        f = np.linalg.solve(btv.T, r - bpos[b])
        if np.abs(f - np.round(f)).max() < 1e-3: return b
cand = {sub(i): i for i in d["at"]}
print("# camadas do volume:", [(s, sorted({round(np.linalg.norm(R), 3) for R in neigh(0) if SH[s][0] < np.linalg.norm(R) < SH[s][1]})) for s in range(3)])
print("# Hb usado por subrede:", cand)
for nsh in (2, 3):
    A, y = [], []
    for b, i in sorted(cand.items()):
        T = tl(d[f"Hb{i}"][3:6, 3:6]); Bs = basis(b, nsh)
        iu = [(0, 0), (1, 1), (0, 1), (0, 2), (1, 2)]
        A += [[Bk[p] for Bk in Bs] for p in iu]; y += [T[p] for p in iu]
    A, y = np.array(A), np.array(y); beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    print(f"\n## {nsh} camadas: β = {np.round(beta, 4)} eV; ‖A_tl‖ volume = {np.linalg.norm(y):.4f}, resíduo ‖y − Aβ‖ = {np.linalg.norm(y - A @ beta):.2e}")
    # conferência N11 (DFT slab, diagnóstico)
    for i in (30, 32, 0, 1):
        mis = d[f"mis{i}"]; yv = None
        cf = -sum((beta[[s for s in range(nsh) if SH[s][0] < np.linalg.norm(m) < SH[s][1]][0]] * U(m) for m in mis), np.zeros((3, 3)))
        dft = tl((d[f"Hs{i}"] - d[f"Hb{i}"])[3:6, 3:6])
        cs = np.sum(cf * dft) / (np.linalg.norm(cf) * np.linalg.norm(dft))
        print(f"  N11 át {i} ({len(mis)} aus, |d| {np.round(np.linalg.norm(mis, axis=1), 2)}): ‖ΔH_cf‖ {np.linalg.norm(cf):.4f} × ‖DFT Hs−Hb‖ {np.linalg.norm(dft):.4f}; cos tensorial {cs:+.2f}")
    # face do L8: ausentes = sítios do volume com n2 ≥ 8 (como a12_termo.py)
    aus = read_sitepot(D + "campo_semiinf/termo/termo_aus.sitepot")[1]
    TEN = {}
    for a in range(13, 25):
        i = a - 1; Ri = posf[i]
        try: yv, _ = bond_normal(i)
        except AssertionError: yv = np.array([0., 1, 0])
        mis = []
        # reconstrói ausentes pela pilha do volume: posf[a] = posb[a%3] + (a//3)·latb[1] (filme = n2 0..7)
        latb1 = posf[3] - posf[0]
        for n2 in range(8, 11):
            for b in range(3):
                for n1 in range(-2, 3):
                    for n3 in range(-2, 3):
                        R = posf[b] + n2 * latb1 + n1 * latf[0] + n3 * latf[2] - Ri
                        if np.linalg.norm(R) < SH[nsh - 1][1]: mis.append(R)
        cf = -sum((beta[[s for s in range(nsh) if SH[s][0] < np.linalg.norm(m) < SH[s][1]][0]] * U(m) for m in mis), np.zeros((3, 3)))
        ta = tl(aus[i][1][3:6, 3:6]) if i in aus else np.zeros((3, 3)); R = cf - ta
        if np.abs(R).max() > 1e-3: TEN[i] = R
        if not mis: continue
        for nm, M in (("ΔH_cf (total)", cf), ("termo_aus VNA+xc", ta), ("R = ΔH_cf − VNA+xc", R)):
            cp = 1.5 * yv @ M @ yv; cs = np.sum(M * U(yv)) / (np.linalg.norm(M) * np.linalg.norm(U(yv)) + 1e-30)
            print(f"  L8 át {a} ({len(mis)} aus, |d| {np.round(sorted(np.linalg.norm(mis, axis=1)), 2)}) {nm:20s}: c_pred = {cp:+.4f} eV; cos com (ŷŷ−I/3) {cs:+.2f}; autoval {np.round(np.linalg.eigvalsh(M), 4)}")
    np.save(D + f"barreira/a20/c_vol/R_cf_{nsh}sh.npy", TEN, allow_pickle=True)
    print(f"  tensores R gravados (átomos {[k + 1 for k in TEN]}): c_vol/R_cf_{nsh}sh.npy")
