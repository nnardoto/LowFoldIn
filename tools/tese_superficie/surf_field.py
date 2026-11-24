#!/usr/bin/env python3
"""A8, sistemas 2 e 3: E_loc na superfície semi-infinita Te(10-10) pelo delta V_H completo.

uso: surf_field.py <dir> <prefixo_pairs> <P> <M> [tag]
  P = pontos por período do volume na normal (a2y/P), M = períodos de volume abaixo da costura.
Pilha na normal (y, cartesiana; x = a, z = c periódicos):
  [M períodos de delta rho do VOLUME (vol/vol_k4.rho_pairs.dat), cópias roladas de a/2]
  [delta rho do sistema (rho_pairs do surfgreen) de y_costura até o vácuo]
y_costura = max y dos átomos não centrais + 7 Bohr (acima dela nenhum par costura-costura falta).
Poisson aberto por G_par (a8lib.open_poisson), G=0 com campo macroscópico nulo no volume.
"""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a8lib as L

# ROOT = pasta com bin/rhomap, vol/ e as pastas dos sistemas (padrão: diretório corrente)
ROOT = os.environ.get("A8_ROOT", os.getcwd())
RHOMAP = os.environ.get("A8_RHOMAP", os.path.join(ROOT, "bin", "rhomap"))
VOLP = os.path.join(ROOT, "vol", os.environ.get("A8_VOL", "vol_k4") + ".rho_pairs.dat")
PAO = L.PAO
NX, NZ = 48, 60
d_sys, pre, P, M = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
tag = sys.argv[5] if len(sys.argv) > 5 else f"P{P}_M{M}"
pairs = os.path.join(ROOT, d_sys, pre + ".rho_pairs.dat")


def lattice(pf):
    lat = {}
    for l in open(pf):
        if l.startswith("lattice "):
            t = l.split(); lat[int(t[1])] = np.array(t[2:5], float)
    return lat


def central_positions(pf):
    ri = set()
    for l in open(pf):
        if l.startswith("pair "):
            t = l.split(); ri.add(tuple(np.round(np.array(t[12:15], float), 8)))
    return np.array(sorted(ri, key=lambda r: r[1]))


def run_rhomap(pf, out, org, a, dy, c, ny):
    nml = os.path.join(ROOT, d_sys, out + ".nml")
    open(nml, "w").write(
        f"&RHOMAP source = 'pairs', pairs_file = '{pf}', pairs_window = 'occ', spin_factor = 2.0, "
        f"output_prefix = '{os.path.join(ROOT, d_sys, out)}' /\n"
        f"&SPECIES nspec = 1, label = 'Te', pao_file = '{PAO}', basis = 's3p3d2f1', rcut_bohr = 7.0 /\n"
        f"&GRID mode = 'cube', units = 'ang', origin = {org[0]:.12f}, {org[1]:.12f}, {org[2]:.12f},\n"
        f"  du = {a/NX:.15f}, 0.0, 0.0, dv = 0.0, {dy:.15f}, 0.0, dw = 0.0, 0.0, {c/NZ:.15f},\n"
        f"  nu = {NX}, nv = {ny}, nw = {NZ} /\n")
    if not os.path.exists(os.path.join(ROOT, d_sys, out + ".cube.bin")):
        with open(os.path.join(ROOT, d_sys, out + ".log"), "w") as lg:
            subprocess.run([RHOMAP, nml], stdout=lg, stderr=subprocess.STDOUT, check=True,
                           env=dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "24")))
    return L.read_cube(os.path.join(ROOT, d_sys, out))


# ---------------- geometria (Å)
lat_s, lat_v = lattice(pairs), lattice(VOLP)
a_s, c_s = lat_s[1][0], lat_s[3][2]
a_v, c_v = lat_v[1][0], lat_v[3][2]
a2v = lat_v[2] / 6.0                                   # TeBulk_x1_L6: 6 colunas
A2Y = a2v[1]
dy = A2Y / P
sites_s = L.read_sites(pairs)[0] * L.B
sites_v = L.read_sites(VOLP)[0] * L.B
cent = central_positions(pairs)
cent_v = central_positions(VOLP)
ncen = np.array([not np.any(np.all(np.abs(s - cent) < 1e-3, axis=1)) for s in sites_s])
# (dobra x, z de cent para comparar com sites dobrados)
cf = cent.copy(); cf[:, 0] = np.mod(cf[:, 0] + 1e-4, a_s) - 1e-4; cf[:, 2] = np.mod(cf[:, 2] + 1e-4, c_s) - 1e-4
ncen = np.array([not np.any(np.all(np.abs(s - cf) < 1e-3, axis=1)) for s in sites_s])
y_seam = sites_s[ncen, 1].max()
y_st = y_seam + 7.0 * L.B + 1e-3
y_top = cent[:, 1].max() + 7.0 * L.B + 4.0
ny_s = int(np.ceil((y_top - y_st) / dy)) + 1
# translação T: sistema = volume + T (casar os 3 átomos não centrais mais altos com sítios do volume)
top3 = sites_s[ncen][np.argsort(sites_s[ncen, 1])[-3:]]
T = None
for rv in sites_v:
    Tc = top3[-1] - rv
    ok = True
    for r in top3:
        d = r - Tc - sites_v
        # módulo a (x), c (z) e a2 (x, y)
        n2 = np.round(d[:, 1] / A2Y)
        d = d - n2[:, None] * a2v
        d[:, 0] -= np.round(d[:, 0] / a_v) * a_v; d[:, 2] -= np.round(d[:, 2] / c_v) * c_v
        if np.min(np.linalg.norm(d, axis=1)) > 2e-3:
            ok = False; break
    if ok:
        T = Tc; break
assert T is not None, "sem translação sistema <-> volume"
# fatia de volume: [y_st - A2Y, y_st] no sistema -> coordenadas do volume (dentro da região completa)
org_v = np.array([0.0, y_st - A2Y, 0.0]) - T
yv_lo = sites_v[:, 1].min()
cv_y = cent_v[:, 1]
lo_ok = cv_y.min() - 1.1 + 7 * L.B  # aproximado; ajusta com n * a2 até caber no meio dos centrais
mid = 0.5 * (cv_y.min() + cv_y.max())
n = int(np.round((mid - (org_v[1] + A2Y / 2)) / A2Y))
org_v = org_v + n * a2v
print(f"[{tag}] a = {a_s:.6f}, c = {c_s:.6f} Å; dy = {dy:.4f} Å; y_costura = {y_st:.3f} Å (não central mais alto {y_seam:.3f});"
      f" T = {np.round(T, 5)}; fatia do volume em y = [{org_v[1]:.3f}, {org_v[1]+A2Y:.3f}] (centrais do volume "
      f"{cv_y.min():.2f}..{cv_y.max():.2f})")

rs, hs = run_rhomap(pairs, f"cube_{pre}_P{P}", [0.0, y_st, 0.0], a_s, dy, c_s, ny_s)
rv, hv = run_rhomap(VOLP, f"cube_{os.path.basename(VOLP).split(".")[0]}_for_{pre}_P{P}", org_v, a_v, dy, c_v, P + 1)
dxs, dzs = hs["du_bohr"][0], hs["dw_bohr"][2]; dyb = dy / L.B
xv_, den = L.atomic_density(); fat = L.make_rho_atom(xv_, den)
o_s = hs["origin_bohr"]; o_v = hv["origin_bohr"]
ra_s = L.rho_atom_grid(fat, sites_s / L.B, o_s[0], o_s[1], o_s[2], dxs, dyb, dzs, NX, ny_s, NZ, a_s / L.B, c_s / L.B)
ra_v = L.rho_atom_grid(fat, sites_v / L.B, o_v[0], o_v[1], o_v[2], hv["du_bohr"][0], dyb, hv["dw_bohr"][2],
                       NX, P + 1, NZ, a_v / L.B, c_v / L.B)
dr_s = rs - ra_s
dr_v = rv - ra_v
dV = dxs * dyb * dzs
jump_rho = np.abs(rs[:, 0, :] - rv[:, P, :]).max()
jump_d = np.abs(dr_s[:, 0, :] - dr_v[:, P, :]).max()
print(f"[{tag}] salto na costura: max|rho_sis - rho_vol| = {jump_rho:.3e} e/Bohr^3 (max rho {rs.max():.2f});"
      f" max|drho_sis - drho_vol| = {jump_d:.3e}")
print(f"[{tag}] carga: período de volume int drho = {dr_v[:, :P].sum()*dV:.3e} e; sistema (y > costura) "
      f"int drho = {dr_s.sum()*dV:.4e} e, int rho = {rs.sum()*dV:.4f}, int rho_at = {ra_s.sum()*dV:.4f}")
slab = dr_v[:, :P]
stack = np.concatenate([np.roll(slab, ((M - 1 - m) * NX) // 2, axis=0) for m in range(M)] + [dr_s], axis=1)
y0 = o_s[1] - M * P * dyb
V, dVy, gx, gz, Q = L.open_poisson(stack, dxs, dyb, dzs)
E0b = dVy[0, :, 0].real
per = [E0b[k * P:(k + 1) * P].mean() for k in range(M)]
cfield = -np.mean(per[1:M])                        # campo macroscópico nulo no volume (exceto o período do fundo)
print(f"[{tag}] G=0: <dV0/dy> por período do volume (V/Å, antes de c): {np.round(np.array(per)*L.EF, 4)}; c = {cfield*L.EF:.4f}")
print(f"[{tag}] campo normal G=0 no topo (vácuo): {(E0b[-1] + cfield)*L.EF:.4e} V/Å; carga líquida da pilha "
      f"{Q[-1]*(dxs*NX)*(dzs*NZ):.3e} e por célula de superfície")
Ef = L.fields_at(V, dVy, gx, gz, o_s[0], y0, o_s[2], dyb, cent / L.B, cfield)
# camadas: colunas de 3 átomos centrais ordenados em y (de baixo para cima)
lay = np.arange(len(cent)) // 3
nl = lay.max() + 1
E0 = np.array([Ef[lay == j][:, [0, 3, 2]].mean(0) for j in range(nl)])
mod = np.array([np.linalg.norm(Ef[lay == j][:, [0, 3, 2]], axis=1).mean() for j in range(nl)])
np.set_printoptions(suppress=True, precision=4, linewidth=160)
print(f"[{tag}] max|E_y(DF6) - E_y(analítico)| = {np.abs(Ef[:,1]-Ef[:,3]).max():.2e} V/Å")
print(f"[{tag}] E0(j) (x, y, z) V/Å e <|E|>, j = camada a partir do topo (1 = externa):")
for j in range(nl - 1, -1, -1):
    d7 = cent[lay == j, 1].min() - y_seam
    print(f"   j={nl-j}  y~{cent[lay==j,1].mean():6.2f} Å  E0 = {E0[j]}  <|E|> = {mod[j]:.4f}"
          + ("   (a < 7 Bohr da costura: não interpretar)" if d7 < 7 * L.B + 1.0 else ""))
np.savez(os.path.join(ROOT, d_sys, f"E_{pre}_{tag}.npz"), T=T, a2v=a2v, av=a_v, cv=c_v, pos=cent, E=Ef[:, [0, 3, 2]], Ey_fd6=Ef[:, 1],
         E0=E0, mod=mod, y_seam=y_seam, y_st=y_st, cfield=cfield * L.EF)
with open(os.path.join(ROOT, d_sys, f"E_{pre}_{tag}.dat"), "w") as f:
    f.write(f"# A8 E_loc (V/Å) = +(1/e) grad dV_H; {pre}, P={P}, M={M}; x = a, y = normal (vácuo em +y), z = c\n")
    f.write("# x y z (Å)  E_x E_y E_z (E_y analítico)  E_y(DF6)  camada_do_topo\n")
    for r, e, ea, l in zip(cent, Ef[:, [0, 3, 2]], Ef[:, 1], lay):
        f.write(f"{r[0]:10.5f} {r[1]:10.5f} {r[2]:10.5f}  {e[0]:9.5f} {e[1]:9.5f} {e[2]:9.5f}  {ea:9.5f}  {nl-l}\n")
