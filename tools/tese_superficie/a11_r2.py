#!/usr/bin/env python3
"""A11, rodada 2 (cubos do run com VNA na malha, fo2/):
1) VNA em espaço real: VNA_cube(r) = Σ_a f(|r − R_a|), f em B-splines cúbicas uniformes (f ≡ 0 para r ≥ R_c, C²),
   mínimos quadrados; R_c cresce até o resíduo rms chegar ao arredondamento do cubo (4 algarismos: ulp/√12 por ponto).
2) δV_xc(r) = V_xc[ρ_vol − Σ_aus ρ_at,k] − V_xc[ρ_vol], PBE (a11_pbe.py, diferenças centrais na malha do OpenMX),
   ρ_vol = tden (fo2), ρ_at = valence.charge.density do Te7.0.pao renormalizada a 16 e (a8lib), supercélula 3×3×3.
   Validação: V_xc[tden] × vxc0.cube.
3) ΔH^aus = −Σ_k⟨μ|f(|r−d_k|)|ν⟩ + ⟨μ|δV_xc|ν⟩ contra ΔH^DFT (a11_blocos.npz); δV_H da rodada 1 só reportado."""
import sys, numpy as np
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
# raízes configuráveis (padrão: layout original da campanha)
LFI_TESE = _os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/")
PAO = _os.environ.get("LFI_PAO", _os.path.join(_HERE, "..", "..", "test", "data", "pao", "Te7.0.pao"))
D = LFI_TESE + "vizinhos/"
sys.path.insert(0, _HERE)
from a11_pbe import GridFD, vxc
import a8lib
HA = 27.211386245988; B2A = 0.529177210903
out = open(D + "a11_r2_saida.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.write(s + "\n")
def readcube(f):
    L = open(f).read().split("\n")
    nat = int(L[2].split()[0]); org = np.array(L[2].split()[1:4], float)
    n = np.array([int(L[3 + k].split()[0]) for k in range(3)])
    ax = np.array([L[3 + k].split()[1:4] for k in range(3)], float)
    tau = np.array([L[6 + a].split()[2:5] for a in range(nat)], float)
    tok = " ".join(L[6 + nat:]).split()
    return np.array(tok, float).reshape(n), org, n, ax, tau, tok
vna, org, n, ax, tau, tok = readcube(D + "fo2/TeBulk.vna.cube")
A = ax * n[:, None]; Om = abs(np.linalg.det(A)); nat = len(tau)
idx = np.stack(np.meshgrid(*[np.arange(k) for k in n], indexing="ij"), -1).reshape(-1, 3)
rpts = org + idx @ ax

# ---------- 1. VNA por B-splines em espaço real ----------
mant = np.array([len(t.split("E")[0].replace("-", "").replace(".", "")) for t in tok[:5]])
ulp = np.array([10.0 ** (int(t.split("E")[1]) - (mant[0] - 1)) for t in tok])
rround = np.sqrt(np.mean(ulp ** 2 / 12))
say(f"# [1] VNA: cubo com {mant[0]} algarismos; arredondamento rms = {rround * HA:.2e} eV (ulp/√12)")
def bspl(x):
    ax_ = np.abs(x)
    return np.where(ax_ < 1, 2 / 3 - ax_ ** 2 + ax_ ** 3 / 2, np.where(ax_ < 2, (2 - ax_) ** 3 / 6, 0.0))
HK = 0.12  # passo dos nós (bohr); com 0,2 o resíduo para em 2,75× o arredondamento (resolução perto do caroço), com 0,12 e 0,08 em ~1,0×
dists = []   # (índice do ponto, distância) para todas as imagens até 12 bohr
for t0 in tau:
    for sh in np.array(np.meshgrid(*[range(-3, 4)] * 3)).reshape(3, -1).T:
        d = np.linalg.norm(rpts - t0 - sh @ A, axis=1); s = d < 12.0
        dists.append((np.nonzero(s)[0], d[s]))
pi_all = np.concatenate([p for p, _ in dists]); d_all = np.concatenate([d for _, d in dists])
def fit(Rc):
    K = int(np.floor((Rc - 2 * HK) / HK + 1e-9)); cen = HK * np.arange(-1, K + 1)
    s = d_all < Rc; P, Dd = pi_all[s], d_all[s]
    M = np.zeros((len(rpts), len(cen)))
    for j, c in enumerate(cen):
        b = bspl((Dd - c) / HK); np.add.at(M[:, j], P, b)
    c, *_ = np.linalg.lstsq(M, vna.ravel(), rcond=1e-10)
    res = vna.ravel() - M @ c
    return cen, c, res
rr = np.linspace(0, 20, 2001)
say("# R_c (bohr)  resíduo rms (eV)  máx (eV)  rms/arredondamento")
best = None
for Rc in np.arange(4.0, 12.01, 0.5):
    cen, c, res = fit(Rc)
    rat = np.sqrt(np.mean(res ** 2)) / rround
    say(f"#   {Rc:5.2f}  {np.sqrt(np.mean(res ** 2)) * HA:.2e}  {np.abs(res).max() * HA:.2e}  {rat:7.2f}")
    if best is None and rat < 1.5:   # critério: rms ≤ 1,5 × arredondamento (nível do arredondamento)
        best = (Rc, cen, c, res)
if best is None: best = (Rc, cen, c, res); say("#   (critério não atingido até 12 bohr; usa o maior R_c)")
Rc, cen, c, res = best
fV = (bspl((rr[:, None] - cen[None, :]) / HK) @ c)
say(f"# R_c escolhido = {Rc:.2f} bohr (1º com resíduo rms ≤ 1,5× o arredondamento)")
g0 = np.load(D + "a11_vat_TeBulk.npz")["vat_VNA_g0"]
say("# f(r) ajuste × v_at g0 (eV): " + " ".join(f"r={x}:{fV[x * 100] * HA:+.3f}/{g0[x * 100] * HA:+.3f}" for x in (1, 2, 3, 4, 5, 6, 7, 8)))
say(f"#   máx|f − g0| em r ∈ [3, 9] bohr: {np.abs(fV - g0)[300:900].max() * HA:.3f} eV; ∫f d³r = {4 * np.pi * np.trapezoid(fV * rr ** 2, rr):.2f} vs Ω⟨VNA⟩/3 = {Om * vna.mean() / nat:.2f} Ha·bohr³")
np.savez(D + "a11_r2_fvna.npz", r=rr, f=fV, Rc=Rc)

# ---------- 2. δV_xc pelo funcional ----------
rho = readcube(D + "fo2/TeBulk.tden.cube")[0]; vx0 = readcube(D + "fo2/TeBulk.vxc0.cube")[0]
dmin = np.full(len(rpts), 1e9)
for t0 in tau:
    for sh in np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).reshape(3, -1).T:
        dmin = np.minimum(dmin, np.linalg.norm(rpts - t0 - sh @ A, axis=1))
g1 = GridFD(A, n)
for gga in (True, False):
    dd = ((vxc(rho, g1, gga) - vx0) * HA).ravel()
    say(f"# [2] V_xc[tden] {'PBE' if gga else 'LDA'} − vxc0 (eV), rms/máx por d_núcleo: " + " ".join(
        f"{a}–{b}:{np.sqrt(np.mean(dd[(dmin >= a) & (dmin < b)] ** 2)):.3f}/{np.abs(dd[(dmin >= a) & (dmin < b)]).max():.2f}"
        for a, b in ((0, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 99))))
xv, den = a8lib.atomic_density(pao=PAO, zval=16.0); frho = a8lib.make_rho_atom(xv, den)
NS = 3; nS = n * NS; AS = A * NS; gS = GridFD(AS, nS)
rhoS = np.tile(rho, (NS, NS, NS))
idxS = np.stack(np.meshgrid(*[np.arange(k) for k in nS], indexing="ij"), -1).reshape(-1, 3)
rS = org + idxS @ ax; dVol = Om / rho.size
iAS = np.linalg.inv(AS)
def mimg(x):   # vetores mínima imagem na supercélula
    f = x @ iAS; return (f - np.round(f)) @ AS
VXC0 = {g: vxc(rhoS, gS, g) for g in (True, False)}

# ---------- PAO ----------
P = open(PAO).read().split("\n")
def block(tag):
    i0 = [k for k, l in enumerate(P) if l.strip().startswith("<" + tag)][0]
    i1 = [k for k, l in enumerate(P) if l.strip().startswith(tag + ">")][0]
    return np.array([l.split() for l in P[i0 + 1:i1]], float)
R0 = block("pseudo.atomic.orbitals.L=0"); R1 = block("pseudo.atomic.orbitals.L=1")
rp = R0[:, 1]; Rs = R0[:, 2:5].T; Rpz = R1[:, 2:5].T; rcut = rp[-1]
def orbitals(X, Y, Z):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2); rs = np.where(r > 1e-12, r, 1e-12)
    ph = [np.interp(r, rp, Rs[z], right=0) / np.sqrt(4 * np.pi) for z in range(3)]
    cc = np.sqrt(3 / (4 * np.pi))
    for z in range(3):
        R = np.interp(r, rp, Rpz[z], right=0); ph += [cc * R * X / rs, cc * R * Y / rs, cc * R * Z / rs]
    return np.array(ph)
def matel_fine(dk, vr, h=0.10, chunk=40):
    x = np.arange(-rcut, rcut + h / 2, h); V = np.zeros((12, 12)); S = np.zeros((12, 12))
    for i0 in range(0, len(x), chunk):
        X, Y, Z = (a.ravel() for a in np.meshgrid(x[i0:i0 + chunk], x, x, indexing="ij"))
        s = X ** 2 + Y ** 2 + Z ** 2 < rcut ** 2; X, Y, Z = X[s], Y[s], Z[s]
        ph = orbitals(X, Y, Z); S += ph @ ph.T * h ** 3
        pot = -sum(np.interp(np.sqrt((X - d[0]) ** 2 + (Y - d[1]) ** 2 + (Z - d[2]) ** 2), rr, vr) for d in dk)
        V += (ph * pot) @ ph.T * h ** 3
    return S, V * HA

dat = np.load(D + "a11_blocos.npz")
SP = list(range(12)); s1p1 = [0, 3, 4, 5]
def metrics(dH, S):
    cc = np.sum(S[np.ix_(s1p1, s1p1)] * dH[np.ix_(s1p1, s1p1)]) / np.sum(S[np.ix_(s1p1, s1p1)] ** 2)
    pb = dH[3:6, 3:6]; ap = pb - np.trace(pb) / 3 * np.eye(3); w, U = np.linalg.eigh(ap)
    k = np.argmax(np.abs(w)); oth = [j for j in range(3) if j != k]
    return cc, w[k] - w[oth].mean(), U[:, k], dH[0, 3:6], np.linalg.norm(dH[np.ix_([0, 1, 2], list(range(3, 12)))])
r1 = np.load(D + "a11_termo_TeBulk.npz")   # δV_H da rodada 1 (só reportado)
say("# métricas (eV): escalar = tr[SΔH]/tr[S²] em {s1,p1}; σ−π = λ_eixo − média(outros) do bloco p1 sem traço; |v| = |ΔH[s1,p1]|; ‖sp‖ Frobenius s(3ζ)–p(3ζ); cos contra a DFT")
say("# át n | fonte            | escalar   σ−π     | |v|     | ‖sp‖   | cos(eixo) cos(v)")
for i in (0, 1, 30, 32, 2, 31):
    dk = dat[f"mis{i}"] / B2A
    dHd = (dat[f"Hs{i}"] - dat[f"Hb{i}"])[np.ix_(SP, SP)]; Sd = dat[f"S{i}"][np.ix_(SP, SP)]
    md = metrics(dHd, Sd)
    rows = [("DFT", dHd)]
    if len(dk):
        Sf, Vna = matel_fine(dk, fV)
        # posição do átomo i no cubo do volume (mesmo referencial cartesiano; translação t do diag_ausentes)
        Ri = (dat["pos"][i] - dat["t"]) / B2A
        rel = mimg(rS - Ri); near = np.sum(rel ** 2, 1) < rcut ** 2
        ph = orbitals(*rel[near].T)
        Sg_ = ph @ ph.T * dVol
        rsub = rhoS.ravel().copy(); dk_core = np.full(len(rS), 1e9)
        for d in dk:
            dv = np.linalg.norm(mimg(rS - Ri - d), axis=1); m = dv < a8lib.RCUT_ATOM
            rsub[m] -= frho(dv[m]); dk_core = np.minimum(dk_core, dv)
        neg = np.sum(rsub < 1e-6)
        rsub = rsub.reshape(nS)
        Vx = {}
        for gga in (True, False):
            dvx = (vxc(rsub, gS, gga) - VXC0[gga]).ravel()[near]
            Vx[gga] = (ph * dvx) @ ph.T * dVol * HA
        # parte sensível ao caroço parcial (pcc, só em < 2 bohr dos núcleos): pontos a < 2 bohr de k ou de i
        dvx = (vxc(rsub, gS, True) - VXC0[True]).ravel()[near]
        core = (dk_core[near] < 2.0) | (np.linalg.norm(rel[near], axis=1) < 2.0)
        Vcore = (ph[:, core] * dvx[core]) @ ph[:, core].T * dVol * HA
        # checagem de quadratura: VNA na malha do cubo × malha fina
        pot = -sum(np.interp(np.linalg.norm(rel[near] - d, axis=1), rr, fV) for d in dk)
        Vna_c = (ph * pot) @ ph.T * dVol * HA
        say(f"#   át {i}: |S_malha − S_scfout| {np.abs(Sg_ - Sd).max():.1e}, |S_fina − S_scfout| {np.abs(Sf - Sd).max():.1e};"
            f" VNA malha×fina |Δ|máx {np.abs(Vna_c - Vna).max():.4f} eV; pontos com ρ_vol − Σρ_at < 1e-6: {neg}")
        rows += [("VNA (ajuste)", Vna), ("δV_xc PBE", Vx[True]), ("δV_xc LDA", Vx[False]), ("  xc <2 bohr núcl.", Vcore),
                 ("δV_H (rodada 1)", r1[f"dH_{i}_dVH_const"]), ("TOTAL VNA+xc", Vna + Vx[True]),
                 ("resto DFT−TOTAL", dHd - Vna - Vx[True])]
    else:
        rows += [("TOTAL (0 ausentes)", np.zeros((12, 12))), ("resto DFT−TOTAL", dHd)]
    for k, M in rows:
        cc, sp_, axv, vec, nsp = metrics(M, Sd)
        ca = abs(axv @ md[2]); cv = vec @ md[3] / (np.linalg.norm(vec) * np.linalg.norm(md[3]) + 1e-30)
        say(f"{i:3d} {len(dk)} | {k:18s} | {cc:+.4f} {sp_:+.4f} | {np.linalg.norm(vec):.4f} | {nsp:.4f} | {ca:.2f} {cv:+.2f}")
out.close()
