#!/usr/bin/env python3
"""Sistema 1 (volume): delta rho do LowFoldIn -> E_loc; periódico x aberto (controle) x DFT (TeBulk.dden.cube)."""
import sys, numpy as np
sys.path.insert(0, sys.path[0]); import a8lib as L
np.set_printoptions(suppress=True, precision=4, linewidth=150)
NPER = int(sys.argv[1]) if len(sys.argv) > 1 else 40
rho, hdr = L.read_cube(f"vol/volslab_h{NPER}")
o, dx, dy, dz = hdr["origin_bohr"], hdr["du_bohr"][0], hdr["dv_bohr"][1], hdr["dw_bohr"][2]
nx, ny, nz = rho.shape; ax, cz = nx * dx, nz * dz
sites, lat = L.read_sites("vol/vol_k4.rho_pairs.dat")
xv, den = L.atomic_density(); fat = L.make_rho_atom(xv, den)
ra = L.rho_atom_grid(fat, sites, o[0], o[1], o[2], dx, dy, dz, nx, ny, nz, ax, cz)
dr = rho - ra
P = NPER
# translação a2 = (-a/2, a2y): n(x, y + a2y) = n(x + a/2, y) -> plano P = plano 0 rolado de -nx/2
j1 = np.abs(rho[:, P, :] - np.roll(rho[:, 0, :], -nx // 2, axis=0)).max()
j2 = np.abs(rho[:, 2 * P, :] - rho[:, 0, :]).max()
print(f"h = a2y/{P} = {dy*L.B:.4f} Å; salto de rho (1 período, 2 períodos): {j1:.2e} {j2:.2e} e/Bohr^3 (max rho {rho.max():.2f})")
dV = dx * dy * dz
for k in range(2):
    print(f"período {k}: int rho = {rho[:, k*P:(k+1)*P].sum()*dV:.6f}, int rho_at = {ra[:, k*P:(k+1)*P].sum()*dV:.6f}, "
          f"int drho = {dr[:, k*P:(k+1)*P].sum()*dV:.2e} e (3 átomos: 48)")
# átomos dentro da caixa de 2 períodos
box = []
for s in sites:
    for m in range(-2, 3):
        for q in range(-2, 3):
            r = s + np.array([m * ax, 0, q * cz])
            if 0 <= r[0] - o[0] < ax and 0 <= r[2] - o[2] < cz and o[1] <= r[1] < o[1] + 2 * P * dy:
                box.append(r)
box = np.unique(np.round(np.array(box), 6), axis=0); box = box[np.argsort(box[:, 1])]
# (1) referência periódica 3D (caixa retangular de 2 períodos)
Ep = L.periodic_fields(dr[:, :2 * P], dx, dy, dz, o[0], o[1], o[2], box)
# (2) máquina aberta: M períodos empilhados (cópia rolada), átomos dos períodos centrais
M = 8
slab = dr[:, :P]
stack = np.concatenate([np.roll(slab, (m * nx) // 2, axis=0) for m in range(M)], axis=1)   # período m: n(x,y+m a2y)=n(x+m a/2,y)
V, dVy, gx, gz, Q = L.open_poisson(stack, dx, dy, dz)
# G = 0: campo macroscópico nulo no volume (média de dV0/dy nos períodos 2..M-3 = 0)
mid = slice(2 * P, (M - 2) * P)
cfield = -np.mean(dVy[0, mid, 0].real)
print(f"G=0: campo uniforme somado (zera a média no volume) = {cfield*L.EF:.3e} V/Å; carga total da pilha {Q[-1]:.2e} e/Bohr^2")
cent = []
for k, r in enumerate(box):
    mper = 3
    rr = r.copy(); rr[1] += mper * P * dy - 0.0      # mesmo átomo no período 3 (+3 a2y; x: -3 a/2)
    rr[0] -= 3 * ax / 2
    cent.append(rr)
Eo = L.fields_at(V, dVy, gx, gz, o[0], o[1], o[2], dy, np.array(cent), cfield)
# o mesmo átomo deslocado de m a2: campo igual; os da caixa periódica estão em períodos 0 e 1
print("E_loc periódico (x y z) V/Å e |E|:")
for e in Ep: print(e, f"{np.linalg.norm(e):.4f}")
print("E_loc aberto (x, y_DF6, z, y_analítico):")
for e in Eo: print(e)
d = Eo[:, [0, 3, 2]] - Ep
print(f"max|E_aberto - E_periódico| = {np.abs(d).max():.2e} V/Å; max|E_y(DF6) - E_y(analítico)| = {np.abs(Eo[:,1]-Eo[:,3]).max():.2e}")
np.save(f"vol/E_vol_h{NPER}.npy", np.c_[box * L.B, Ep, Eo])
mod = np.linalg.norm(Ep, axis=1)
print(f"<|E_loc|> = {mod.mean():.4f} V/Å; soma por célula (3 átomos) = {Ep[:3].sum(0)}")
