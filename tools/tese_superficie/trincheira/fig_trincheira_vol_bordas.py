#!/usr/bin/env python3
"""A25: densidade integrada ao longo de c nas janelas [E_F − 0,5; E_F] e [E_F; E_F + 0,5] eV na trincheira.
Linha de cima: DFT (OpenMX, slab com uma cadeia removida, loneChainMapa_k6; mesma leitura do te_trincheira_bordas.py da tese).
Linha de baixo: LowFoldIn só do volume (trench/ do A25: one-shot SCC com v_init = perfil convergido, rho_out na janela,
rhomap 'integrate'; E_F do volume). Mesmo recorte (±17 Å em x a partir da cadeia removida; −10 a +4 Å da normal a partir do
átomo mais externo) e mesma orientação (visto de +c). Escala própria por painel, 10⁻² e/Å².
uso: ~/.venvs/tese-figuras/bin/python fig_trincheira_vol_bordas.py"""
import os, sys, json, numpy as np
from pathlib import Path
FS = Path(os.environ.get("TESE_FIGSCRIPTS", "/home/nathanael/Documentos/Doutorado/Tese-UFES/figuras/scripts")); sys.path.insert(0, str(FS))
import te_cube as tc
import te_trincheira_bordas as tb
TV = Path(os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese")) / "trincheira_vol"
DADOS = FS.parent.parent / "materiais-locais" / "loneChainMapa_k6"
HW, DEP, VAC, h = 17.0, 10.0, 4.0, 0.04
plt = tc.setup_style(usetex=False)
from matplotlib.patches import Circle
xs = np.arange(-HW, HW + h / 2, h); ys = np.arange(-DEP, VAC + h / 2, h)
# ---- DFT ----
cubes = {"val": tc.read_cube(str(DADOS / "bias_-0.50eV" / "loneChain.pden.cube")),
         "cond": tc.read_cube(str(DADOS / "bias_+0.50eV" / "loneChain.pden.cube"))}
info = tc.trench_info(cubes["val"]); xc, yc, ytop, alat = info["center"][0], info["center"][1], info["y_top"], info["alat"]
x0 = float(tc.a_of_x(cubes["val"], xc, tc.j_of_y(cubes["val"], ytop + 1.0), info["period"]))
X, Y = np.meshgrid(xs + xc, ys + ytop)
mD = {k: tc.resample_ij(tc.project_c(c), c, X, Y) for k, c in cubes.items()}
c0 = cubes["val"]; xyz = c0.xyz_ang[c0.Z > 0]
atD = (tc.wrap_x(c0, xyz[:, 0], xc) - xc, xyz[:, 1] - ytop); remD = yc - ytop
# ---- LowFoldIn só volume ----
g = json.load(open(TV / "bordas/geom.json"))
XL, YL = np.meshgrid(xs + g["xc"], ys + g["ytop"])
mL = {k: tb.ler_lfi(str(TV / f"bordas/tr_{w}.map.dat"), XL, YL) for k, w in (("val", "occ"), ("cond", "emp"))}
L = open(TV / "hs/TeFilm_x11_L8_tr.structure.nml").read().split("\n")
lat1 = float([l for l in L if "LatticeVectors(1" in l][0].split("=")[1].split()[0])
j = [k for k, l in enumerate(L) if l.strip() == "&ATOMS"][0]
pos = np.array([[float(x) for x in L[j + 2 + a].split()[1:4]] for a in range(int(L[j + 1]))])
ax_ = (pos[:, 0] - g["xc"] + lat1 / 2) % lat1 - lat1 / 2
atL = (ax_, pos[:, 1] - g["ytop"]); remL = g["yc"] - g["ytop"]
# ---- figura ----
P = {"val": (tc.LARANJA, "tc_val", r"ocupados, $[E_F-0{,}5;\,E_F]$ eV"), "cond": (tc.AZUL, "tc_cond", r"desocupados, $[E_F;\,E_F+0{,}5]$ eV")}
fig, axs = plt.subplots(2, 2, figsize=(7.0, 3.6), constrained_layout=True, sharex=True, sharey=True)
for r, (maps, at, rem, lab) in enumerate([(mD, atD, remD, "DFT (OpenMX, supercélula)"), (mL, atL, remL, "LowFoldIn, só volume")]):
    for cidx, key in enumerate(("val", "cond")):
        ax = axs[r, cidx]; cor, nome, tit = P[key]; m = maps[key]
        im = ax.imshow(m, origin="lower", extent=(xs[0] + x0, xs[-1] + x0, ys[0], ys[-1]), cmap=tc.seq_cmap(cor, nome),
                       vmin=0, vmax=float(m.max()), interpolation="bilinear")
        s = (np.abs(at[0]) < HW + 1) & (at[1] > -DEP - 1) & (at[1] < VAC)
        for x, y in zip(at[0][s], at[1][s]):
            ax.add_patch(Circle((x + x0, y), 0.42, facecolor="0.80", edgecolor="0.35", lw=0.35, zorder=3))
        ax.add_patch(Circle((x0, rem), 1.95, fill=False, ls=(0, (2.5, 1.6)), lw=0.7, edgecolor="0.15", zorder=4))
        for sg, txt in ((-1, "parede esq."), (1, "parede dir.")):
            ax.text(x0 + sg * alat, VAC - 0.4, txt, ha="center", va="top", fontsize=6)
        ax.set_xlim(xs[0] + x0, xs[-1] + x0); ax.set_ylim(ys[0], ys[-1]); ax.set_aspect("equal")
        if r == 0: ax.set_title(tit, fontsize=8)
        if r == 1: ax.set_xlabel(r"posição $x$ (Å), ao longo de $\mathbf{a}_1$")
        if cidx == 0: ax.set_ylabel(f"{lab}\nnormal (Å)", fontsize=8)
        cb = fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.8, pad=0.01, format=tb._fmt_scaled(100))
        cb.set_label(r"$10^{-2}$ e/Å$^2$", fontsize=7); cb.ax.tick_params(labelsize=6)
        # peso por lado: integral do mapa em x < x0 − a/2 (esq.) e x > x0 + a/2 (dir.), camada externa (normal > −4 Å)
        sel = ys > -4.0; mm = m[sel]
        wl = mm[:, xs < -alat / 2].sum(); wr = mm[:, xs > alat / 2].sum()
        k = np.unravel_index(np.argmax(m), m.shape)
        print(f"{lab:28s} {key:4s}: máx {m.max():.4f} e/Å² em x−x_trinch = {xs[k[1]]:+.2f} Å, normal {ys[k[0]]:+.2f} Å;"
              f" peso esq./dir. (|x| > a/2, normal > −4 Å) = {wl / (wl + wr):.2f}/{wr / (wl + wr):.2f}")
fig.savefig(TV / "fig_trincheira_vol_bordas.png", dpi=200, facecolor="white")
print("figura:", TV / "fig_trincheira_vol_bordas.png")
