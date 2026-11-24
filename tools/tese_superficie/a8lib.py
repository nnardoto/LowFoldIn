"""A8: densidade do LowFoldIn -> delta V_H -> E_loc nos núcleos (Te(10-10) semi-infinito).

Unidades internas: Bohr, Hartree, e/Bohr^3. Campo de saída em V/Å, E_loc = +(1/e) grad dV_H
(dV_H = energia do elétron, OpenMX: Poisson de rho - sum rho_atomo).
Malha cartesiana: x (periódica, a), y (normal, aberta), z (periódica, c). Arrays (nx, ny, nz).
"""
import numpy as np
from math import factorial, comb

B = 0.529177210903          # Å/Bohr (CODATA 2018; o OpenMX usa 0.529177249)
HA = 27.211386245988
EF = HA / B                 # Ha/Bohr -> V/Å
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
# raízes configuráveis (padrão: layout original da campanha)
LFI_TESE = _os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/")
PAO = _os.environ.get("LFI_PAO", _os.path.join(_HERE, "..", "..", "test", "data", "pao", "Te7.0.pao"))
RCUT_ATOM = 7.0             # Bohr, Spe_Atom_Cut1 do Te7.0 (GridListAtom do Set_Aden_Grid.c)
ZVAL = 16.0                 # Spe_Core_Charge do Te_PBE19.vps (valence.electron 16: 4d10 5s2 5p4)


# ---------------------------------------------------------------- densidade atômica (OpenMX)
def atomic_density(pao=PAO, zval=ZVAL):
    """<valence.charge.density do .pao, renormalizada como SetPara_DFT.c (soma retangular em log r)."""
    lines = open(pao).read().splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.strip().startswith("<valence.charge.density"))
    rows = []
    for l in lines[i0 + 1:]:
        if l.strip().startswith("valence.charge.density>"):
            break
        rows.append([float(t) for t in l.split()])
    a = np.array(rows)
    xv, rv, den = a[:, 0], a[:, 1], a[:, 2]
    dx = xv[1] - xv[0]
    s = 4 * np.pi * dx * np.sum(den * np.exp(3 * xv))
    return xv, den * zval / s


def make_rho_atom(xv, den):
    h = xv[1] - xv[0]          # malha uniforme em x = log r; Lagrange cúbico de 4 nós

    def f(r):
        r = np.maximum(r, np.exp(xv[0]))
        t = (np.log(r) - xv[0]) / h
        j = np.clip(np.floor(t).astype(int) - 1, 0, len(xv) - 4)
        u = t - j
        out = (-(u - 1) * (u - 2) * (u - 3) / 6 * den[j] + u * (u - 2) * (u - 3) / 2 * den[j + 1]
               - u * (u - 1) * (u - 3) / 2 * den[j + 2] + u * (u - 1) * (u - 2) / 6 * den[j + 3])
        out[r > RCUT_ATOM] = 0.0
        return out
    return f


def rho_atom_grid(f, atoms, x0, y0, z0, dx, dy, dz, nx, ny, nz, ax, cz):
    """sum_a rho_at(|r - R_a|) na malha, com imagens em x (ax) e z (cz). atoms (N,3) Bohr."""
    out = np.zeros((nx, ny, nz))
    xs = x0 + dx * np.arange(nx); ys = y0 + dy * np.arange(ny); zs = z0 + dz * np.arange(nz)
    for R in atoms:
        if R[1] < ys[0] - RCUT_ATOM or R[1] > ys[-1] + RCUT_ATOM:
            continue
        jy = np.nonzero(np.abs(ys - R[1]) <= RCUT_ATOM)[0]
        DY = (ys[jy] - R[1])[None, :, None]
        for mx in range(-3, 4):
            DX = (xs - R[0] - mx * ax)
            if np.abs(DX).min() > RCUT_ATOM:
                continue
            for mz in range(-3, 4):
                DZ = (zs - R[2] - mz * cz)
                if np.abs(DZ).min() > RCUT_ATOM:
                    continue
                r = np.sqrt(DX[:, None, None] ** 2 + DY ** 2 + DZ[None, None, :] ** 2)
                out[:, jy, :] += f(r.ravel()).reshape(r.shape)
    return out


# ---------------------------------------------------------------- E/S do LowFoldIn
def read_cube(prefix):
    hdr = {}
    for l in open(prefix + ".cube.hdr"):
        if l.startswith("#"):
            continue
        t = l.split()
        hdr[t[0]] = np.array(t[1:], float)
    n = hdr["n"].astype(int)
    d = np.fromfile(prefix + ".cube.bin", dtype="<f8").reshape(n[::-1]).transpose(2, 1, 0)
    return d, hdr


def read_sites(pairs_file):
    """Posições (Bohr) dos átomos i e j de rho_pairs.dat (únicas) e vetores de rede (Bohr, flag k)."""
    pos, lat = [], {}
    for l in open(pairs_file):
        if l.startswith("lattice "):
            t = l.split(); lat[int(t[1])] = (np.array(t[2:5], float) / B, int(t[5]))
        elif l.startswith("pair "):
            t = l.split()
            d = np.array(t[8:11], float); ri = np.array(t[12:15], float)
            pos.append(ri); pos.append(ri + d)
    p = np.array(pos)
    for k in (1, 3):                    # eixos periódicos (x = a, z = c, vetores puros): dobra
        v, per = lat[k]
        if per:
            j = int(np.argmax(np.abs(v)))
            L_ = abs(v[j]) * B
            p[:, j] = np.mod(p[:, j] + 1e-4, L_) - 1e-4
    _, idx = np.unique(np.round(p, 4), axis=0, return_index=True)
    p = p[np.sort(idx)] / B             # posições originais (sem arredondar)
    return p, lat


# ---------------------------------------------------------------- Poisson aberto em y (por G_par)
_C = None


def _lagrange_coefs():
    """c[k, n]: ell_k(s) = sum_n c[k,n] s^n para nós s = -1, 0, 1, 2 (intervalo [0,1])."""
    nodes = np.array([-1.0, 0.0, 1.0, 2.0])
    V = np.vander(nodes, 4, increasing=True)          # V[k, n] = s_k^n
    return np.linalg.inv(V).T                          # p(s) = sum_k rho_k ell_k(s)


def _moments(a):
    """m_n(a) = int_0^1 e^{-a(1-s)} s^n ds, n = 0..3 (a >= 0, vetorizado, estável)."""
    a = np.asarray(a, float)
    m = np.zeros((4,) + a.shape)
    small = a < 2.0
    asx = a[small]
    for n in range(4):
        acc = np.zeros_like(asx); term = np.ones_like(asx)
        for k in range(40):
            acc += term * factorial(n) / factorial(n + k + 1)
            term = term * (-asx)
        m[n][small] = acc
    al = a[~small]
    if al.size:
        mm = (1 - np.exp(-al)) / al
        m[0][~small] = mm
        for n in range(1, 4):
            mm = (1 - n * mm) / al
            m[n][~small] = mm
    return m


def open_poisson(n, dx, dy, dz):
    """n(x,y,z) (e/Bohr^3) -> (V_G(y), dV_G/dy(y)) por G_par, kernel aberto (2 pi/G) e^{-G|y-y'|}.

    Integração por produto: n_G(y') cúbico por partes (Lagrange de 4 nós), exponencial exata.
    G = 0: V0 = -2 pi int |y-y'| n0 dy' (sem a parcela linear c y, somada por quem chama),
    dV0/dy = 2 pi (R - L). Retorna também G (|G|) e as frequências.
    """
    nx, ny, nz = n.shape
    F = np.fft.fft2(n, axes=(0, 2)) / (nx * nz)          # (nx, ny, nz)
    gx = 2 * np.pi * np.fft.fftfreq(nx, d=dx); gz = 2 * np.pi * np.fft.fftfreq(nz, d=dz)
    G = np.sqrt(gx[:, None] ** 2 + gz[None, :] ** 2)     # (nx, nz)
    c = _lagrange_coefs()
    a = G * dy
    m = _moments(a)                                       # (4, nx, nz)
    binm = np.array([[comb(n_, j) * (-1) ** j if j <= n_ else 0 for j in range(4)] for n_ in range(4)], float)
    p = np.einsum("nj,jxz->nxz", binm, m)                 # int_0^1 e^{-a s} s^n ds
    wL = dy * np.einsum("kn,nxz->kxz", c, m)              # pesos dos 4 nós, recursão de baixo
    wR = dy * np.einsum("kn,nxz->kxz", c, p)              # recursão de cima
    e = np.exp(-a)
    L = np.zeros_like(F); R = np.zeros_like(F)

    def stencil(i):            # nós do intervalo [i-1, i] (clamp nas bordas)
        j0 = min(max(i - 2, 0), ny - 4)
        return j0, i - 1 - j0  # posição do nó i-1 no estêncil (normalmente 1)
    W = []
    for i in range(1, ny):
        j0, off = stencil(i)
        if off == 1:
            W.append((j0, wL, wR))
        else:  # borda: reconstrói pesos para nós deslocados (s do nó i-1 = 0 sempre)
            nodes = np.arange(4) - off
            Vm = np.vander(nodes.astype(float), 4, increasing=True)
            cc = np.linalg.inv(Vm).T
            W.append((j0, dy * np.einsum("kn,nxz->kxz", cc, m), dy * np.einsum("kn,nxz->kxz", cc, p)))
    for i in range(1, ny):
        j0, wl, _ = W[i - 1]
        L[:, i, :] = e * L[:, i - 1, :] + np.einsum("kxz,xkz->xz", wl, F[:, j0:j0 + 4, :])
    for i in range(ny - 2, -1, -1):
        j0, _, wr = W[i]
        R[:, i, :] = e * R[:, i + 1, :] + np.einsum("kxz,xkz->xz", wr, F[:, j0:j0 + 4, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        pref = np.where(G > 0, 2 * np.pi / G, 0.0)
    V = pref[:, None, :] * (L + R)
    dV = 2 * np.pi * (R - L)
    # G = 0: potencial -2 pi int |y-y'| n0 (primeiros momentos com os mesmos pesos cúbicos, a = 0)
    n0 = np.real(F[0, :, 0])
    ys = dy * np.arange(ny)
    Q = np.zeros(ny); P = np.zeros(ny)
    m0 = np.array([1.0, 1 / 2, 1 / 3, 1 / 4]); m1 = np.array([1 / 2, 1 / 3, 1 / 4, 1 / 5])
    for i in range(1, ny):
        j0 = min(max(i - 2, 0), ny - 4); off = i - 1 - j0
        nodes = np.arange(4) - off
        cc = np.linalg.inv(np.vander(nodes.astype(float), 4, increasing=True)).T
        w0 = dy * cc @ m0; w1 = dy * cc @ m1
        Q[i] = Q[i - 1] + w0 @ n0[j0:j0 + 4]
        P[i] = P[i - 1] + ys[i - 1] * (w0 @ n0[j0:j0 + 4]) + dy * (w1 @ n0[j0:j0 + 4])
    Qt, Pt = Q[-1], P[-1]
    V0 = -2 * np.pi * (ys * Q - P + (Pt - P) - ys * (Qt - Q))
    V[0, :, 0] = V0
    dV[0, :, 0] = 2 * np.pi * ((Qt - Q) - Q)
    return V, dV, gx, gz, Q


def fd6(f, h, axis):
    """Derivada central de 6.ª ordem (pontos interiores; bordas = nan)."""
    c = np.array([-1 / 60, 3 / 20, -3 / 4, 0, 3 / 4, -3 / 20, 1 / 60]) / h
    f = np.moveaxis(f, axis, 0)
    out = np.full(f.shape, np.nan, dtype=f.dtype)
    out[3:-3] = sum(c[k] * f[k:f.shape[0] - 6 + k] for k in range(7))
    return np.moveaxis(out, 0, axis)


def lagr_weights(t, k=6):
    """Pesos de Lagrange (k nós consecutivos) para interpolar em t (em unidades de passo, origem 0)."""
    j0 = int(np.floor(t)) - (k // 2 - 1)
    nodes = j0 + np.arange(k)
    w = np.ones(k)
    for i in range(k):
        for j in range(k):
            if j != i:
                w[i] *= (t - nodes[j]) / (nodes[i] - nodes[j])
    return j0, w


def fields_at(V, dV, gx, gz, x0, y0, z0, dy, atoms, c_field=0.0):
    """E_loc (V/Å) nos átomos: espectral em (x, z), Lagrange 6 pontos em y.
    Rotas: (b) dV_G/dy analítico; (a) DF de 6.ª ordem de V na normal. c_field (Ha/Bohr) = campo uniforme G=0."""
    Vy_a = fd6(V, dy, 1)
    res = []
    for R in atoms:
        ex = np.exp(1j * gx * (R[0] - x0)); ez = np.exp(1j * gz * (R[2] - z0))
        t = (R[1] - y0) / dy
        j0, w = lagr_weights(t)
        sl = slice(j0, j0 + 6)
        def syn(A):
            return np.real(np.einsum("xjz,x,z->j", A[:, sl, :], ex, ez)) @ w
        Ex = syn(1j * gx[:, None, None] * V); Ez = syn(1j * gz[None, None, :] * V)
        Eyb = syn(dV) + c_field; Eya = syn(Vy_a) + c_field
        res.append([Ex, Eya, Ez, Eyb])
    return np.array(res) * EF


def periodic_fields(n, dx, dy, dz, x0, y0, z0, atoms):
    """Referência totalmente periódica (caixa retangular): V_G = 4 pi n_G / G^2, gradiente espectral 3D."""
    nx, ny, nz = n.shape
    F = np.fft.fftn(n) / n.size
    g = [2 * np.pi * np.fft.fftfreq(m, d=h) for m, h in zip(n.shape, (dx, dy, dz))]
    G2 = g[0][:, None, None] ** 2 + g[1][None, :, None] ** 2 + g[2][None, None, :] ** 2
    G2[0, 0, 0] = 1.0
    Vg = 4 * np.pi * F / G2; Vg[0, 0, 0] = 0.0
    out = []
    for R in atoms:
        e = [np.exp(1j * g[k] * (R[k] - (x0, y0, z0)[k])) for k in range(3)]
        ph = np.einsum("i,j,k->ijk", *e)
        out.append([np.real(np.sum(1j * g[k].reshape([-1 if q == k else 1 for q in range(3)]) * Vg * ph))
                    for k in range(3)])
    return np.array(out) * EF
