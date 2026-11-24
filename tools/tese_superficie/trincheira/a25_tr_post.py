#!/usr/bin/env python3
"""A25: pós-processamento da trincheira só-volume (trench/): (1) mapas A(k_c, E) médios nos átomos a < 4 Å do topo de cada coluna
(longe t1 = 2, parede esq. t1 = 4, fundo t1 = 5, parede dir. t1 = 6; coluna = pilha de hélices t1 do hsiocell) lado a lado com a
trincheira LowFoldIn-sobre-DFT-supercélula (trench_slices/, seleção do fig_trench_prelim.py: colunas 3, 5, 6, 7);
(2) carga induzida Mulliken por hélice (dq = s q_DS − q_ref, elétrons por período c) contra a DFT (trench/carga_por_helice.csv,
mesma convenção: + = excesso de elétrons; coluna DFT = t1 + 1, camada DFT = t2 + 4).
uso: a25_tr_post.py [tag]   (lê trench/<tag>.scc_profile.dat e trench/maps/a*_c*.map.dat)"""
import os, re, sys, numpy as np
sys.path.insert(0, os.environ.get("TESE_FIGSCRIPTS", "/home/nathanael/Documentos/Doutorado/Tese-UFES/figuras/scripts"))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import te_cube as tc
R = os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/").rstrip("/") + "/"; TV = R + "trincheira_vol/"; TR = TV + "trench/"
tag = sys.argv[1] if len(sys.argv) > 1 else "tr"
def readstruct(f):
    L = open(f).read().split("\n")
    lat = np.array([[float(x) for x in l.split("=")[1].split()] for l in L if "LatticeVectors" in l])
    i = [k for k, l in enumerate(L) if l.strip() == "&ATOMS"][0]; nat = int(L[i + 1])
    return lat, np.array([[float(x) for x in L[i + 2 + a].split()[1:4]] for a in range(nat)])
latb, posb = readstruct(TV + "hs/TeBulk.structure.nml"); _, pos = readstruct(TV + "hs/TeFilm_x11_L8_tr.structure.nml")
lab = []
for p in pos:
    for b in range(3):
        f = np.linalg.solve(latb.T, p - posb[b])
        if np.abs(f - np.round(f)).max() < 1e-6: lab.append([b, *np.round(f).astype(int)]); break
lab = np.array(lab)
COLS = [(2, "longe (3 col. à esq.)"), (4, "parede esquerda"), (5, "fundo"), (6, "parede direita")]
def selection():
    cen = [int(s) for s in re.findall(r"supercell atom (\d+)", open(TR + f"{tag}_central.nml").read())]
    out = {}
    for c, _ in COLS:
        m = [s for s in cen if lab[s - 1, 1] == c]; ytop = max(pos[s - 1, 1] for s in m)
        out[c] = [s for s in m if pos[s - 1, 1] > ytop - 4.0]
    return out
def loadmap(f):
    d = np.loadtxt(f); k = np.unique(d[:, 0]); e = np.unique(d[:, 1]); return k, e, d[:, 2].reshape(len(k), len(e))
if __name__ == "__main__":
    sel = selection()
    print("# seleção (supercell atom, camada t2):", {c: [(s, int(lab[s - 1, 2])) for s in v] for c, v in sel.items()})
    # ---- carga ----
    dft = np.loadtxt(R + "../trench/carga_por_helice.csv", delimiter=",", skiprows=1)
    dq = {}
    for l in open(TR + f"{tag}.scc_profile.dat"):
        if l.startswith("#"): continue
        f = l.split(); s = int(f[1]); key = (int(lab[s - 1, 1]), int(lab[s - 1, 2])); dq[key] = dq.get(key, 0.0) + float(f[9])
    print("# carga induzida por hélice (elétrons por período c; + = excesso de elétrons): t1 t2 | LowFoldIn-volume | DFT supercélula (col t1+1, camada t2+4)")
    rows = []
    for t2 in (7, 6, 5, 4):
        for t1 in range(11):
            if (t1, t2) not in dq: continue
            dd = dft[(dft[:, 0] == t1 + 1) & (dft[:, 1] == t2 + 4)]; rows.append((t1, t2, dq[(t1, t2)], dd[0, 4] if len(dd) else np.nan))
            print(f"{t1:2d} {t2} | {dq[(t1, t2)]:+.4f} | {rows[-1][3]:+.4f}")
    np.savetxt(TV + f"carga_{tag}.txt", np.array(rows), fmt="%d %d %+.5f %+.5f", header="t1 t2 dq_LowFoldIn_vol dq_DFT (e/c)")
    print(f"# soma central LowFoldIn {sum(dq.values()):+.4f} e/c")
    # ---- mapas ----
    M = {}
    for c, v in sel.items():
        A = [loadmap(TR + f"maps/a{s}_c{lab[s - 1, 2]}.map.dat") for s in v]; k, e = A[0][0], A[0][1]; M[c] = sum(a[2] for a in A) / len(A)
    D = R + "trench_slices/"; txt = open(D + "trenchx_c3_9_central.nml").read()
    at = re.findall(r"^(Te|E)\s+(\S+)\s+(\S+)\s+(\S+)\s+! column (\d+), supercell atom (\d+)", txt, re.M)
    sp = np.array([a[0] for a in at]); y = np.array([float(a[2]) for a in at]); col = np.array([int(a[4]) for a in at]); sid = np.array([int(a[5]) for a in at])
    MD = {}
    for c in (3, 5, 6, 7):
        m = (col == c) & (sp == "Te"); yt = y[m].max(); ii = np.where(m & (y > yt - 4.0))[0]
        A = [loadmap(D + f"a{sid[i]}_c{col[i]}.map.dat") for i in ii]; kd, ed = A[0][0], A[0][1]; MD[c] = sum(a[2] for a in A) / len(A)
    CM = tc.seq_cmap(tc.INDIGO if hasattr(tc, "INDIGO") else "#3f3a8c", "ind")
    vmax = max(max(M[c].max() for c in M), max(MD[c].max() for c in MD)) * 0.5
    fig, axs = plt.subplots(2, 4, figsize=(11, 6.4), sharey=True, sharex=True, constrained_layout=True)
    for j, ((c, t), cd) in enumerate(zip(COLS, (3, 5, 6, 7))):
        for r, (MM, kk, ee, cc, lb) in enumerate([(M[c], k, e, c, "LowFoldIn, só volume"), (MD[cd], kd, ed, cd, "LowFoldIn sobre DFT supercélula")]):
            ax = axs[r, j]; kz = np.linspace(0, 0.5, len(kk))
            im = ax.pcolormesh(kz, ee, MM.T, cmap=CM, vmin=0, vmax=vmax, shading="nearest")
            ax.axhline(0, color="0.4", lw=0.6, ls="--"); ax.set_ylim(-1, 1)
            ax.set_title(f"{t}" if r == 0 else "", fontsize=9)
            if r == 1: ax.set_xlabel(r"$k_c$ ($2\pi/c$)")
            if j == 0: ax.set_ylabel(f"{lb}\n" + r"$E-E_F$ (eV)", fontsize=8)
    fig.colorbar(im, ax=axs, label="A (estados/eV/átomo)", shrink=0.8)
    fig.savefig(TV + f"fig_trincheira_vol_{tag}.png", dpi=150); print("figura:", TV + f"fig_trincheira_vol_{tag}.png")
    np.savez(TV + f"mapas_{tag}.npz", k=k, e=e, **{f"c{c}": M[c] for c in M})
