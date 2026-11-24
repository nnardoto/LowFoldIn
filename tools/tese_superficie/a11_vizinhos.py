#!/usr/bin/env python3
"""A11 — termo de vizinhos ausentes a partir do potencial local do VOLUME do OpenMX, sem parâmetro ajustado.
Por PARTE do potencial (cada uma deconvoluída à parte; o termo é linear, total = soma):
  VNA  = TeBulk.vna.cube (só no run com VNA na malha, level.of.fileout 2)
  dVH  = TeBulk.vhart.cube (δV_H)
  xc   = TeBulk.vxc0.cube, se existir; senão v0 − vhart (− vna, se o v0 contiver o VNA; ver checagem abaixo)
No run original (scf.ProExpn.VNA ligado) o v0.cube = δV_H + V_xc, sem VNA (conferido: v0 − vhart ≈ V_xc, corr. 0,99 com o x LDA).
1) deconvolução v(G) = ṽ_at(|G|) S(G)/Ω por casca |G| (mínimos quadrados Σ S*V/Σ|S|²), v_at(r) pela transformada radial;
2) ΔH^aus_μν(i) = −Σ_k ⟨μ_i| v_at(|r − d_k|) |ν_i⟩ (on-site, s e p com 3 ζ, PAO Te7.0.pao), integração 3D;
3) comparação com ΔH^DFT = H(N11) − H(volume) (a11_blocos.npz, de a11_extrai.py no scicom).
Uso: python a11_vizinhos.py [DIR_DOS_CUBOS] [PREFIXO]   (padrão: esta pasta, TeBulk). Unidades internas bohr/Hartree; saída eV."""
import sys, os, numpy as np
HA = 27.211386245988; B2A = 0.529177210903
D = os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/") + "vizinhos/"
CD = sys.argv[1] if len(sys.argv) > 1 else D
PFX = sys.argv[2] if len(sys.argv) > 2 else "TeBulk"
PAO = os.environ.get("LFI_PAO", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "test", "data", "pao", "Te7.0.pao"))
out = open(D + f"a11_saida_{PFX}.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.write(s + "\n")

# ---------- 0. cubos ----------
def readcube(f):
    L = open(f).read().split("\n")
    nat = int(L[2].split()[0]); org = np.array(L[2].split()[1:4], float)
    n = np.array([int(L[3 + k].split()[0]) for k in range(3)])
    ax = np.array([L[3 + k].split()[1:4] for k in range(3)], float)
    tau = np.array([L[6 + a].split()[2:5] for a in range(nat)], float)
    return np.array(" ".join(L[6 + nat:]).split(), float).reshape(n), org, n, ax, tau
cf = lambda s: os.path.join(CD, f"{PFX}.{s}.cube")
v0, org, n, ax, tau = readcube(cf("v0")); nat = len(tau)
cub = {"v0": v0, "vhart": readcube(cf("vhart"))[0]}
for s in ("vna", "vxc0"):
    if os.path.exists(cf(s)): cub[s] = readcube(cf(s))[0]
parts = {}
if "vna" in cub: parts["VNA"] = cub["vna"]
parts["dVH"] = cub["vhart"]
if "vxc0" in cub:
    parts["xc"] = cub["vxc0"]
    for lab, tot in (("δV_H+V_xc", cub["vhart"] + cub["vxc0"]), ("VNA+δV_H+V_xc", cub["vhart"] + cub["vxc0"] + cub.get("vna", 0))):
        say(f"# checagem: v0 − ({lab}) rms {np.sqrt(np.mean((v0 - tot) ** 2)) * HA:.4f} eV")
else:
    parts["xc"] = v0 - cub["vhart"] - (cub["vna"] if "vna" in cub else 0)   # supõe v0 completo se houver vna.cube
say(f"# partes: {list(parts)}; médias (eV): " + ", ".join(f"{k} {p.mean() * HA:+.3f}" for k, p in parts.items()))

A = ax * n[:, None]; Om = abs(np.linalg.det(A)); Bm = 2 * np.pi * np.linalg.inv(A).T   # linhas = b_i
N = v0.size
m = [np.fft.fftfreq(n[k], 1.0 / n[k]) for k in range(3)]
G = np.stack(np.meshgrid(*m, indexing="ij"), -1) @ Bm; Gn = np.linalg.norm(G, axis=-1)
Sg = np.exp(-1j * ((tau - org) @ G.reshape(-1, 3).T)).sum(0).reshape(n)
Gsph = min(np.pi * n[k] / np.linalg.norm(A[k]) for k in range(3))   # esfera inscrita na caixa |m_k| ≤ n_k/2 (G·a_k = 2π m_k)
key = np.round(Gn, 5); ins = (Gn > 1e-8) & (Gn < Gsph); shells = np.unique(key[ins])
CRIT = 0.05   # descarta casca se ⟨|S|²⟩ < CRIT·máx|S|² (= CRIT·nat²): ṽ não determinado onde S ≈ 0
s2sh = np.array([np.sum(np.abs(Sg[ins & (key == g)]) ** 2) for g in shells])
ngsh = np.array([np.sum(ins & (key == g)) for g in shells])
keep = s2sh / ngsh > CRIT * nat ** 2
say(f"# cubo {n}, Ω = {Om:.2f} bohr³, G_esf = {Gsph:.3f} bohr⁻¹; cascas {len(shells)}, descartadas {np.sum(~keep)}"
    f" (inclui |G| = {', '.join(f'{g:.2f}' for g in shells[~keep][:3])}...); 1ª casca útil {shells[keep][0]:.3f} bohr⁻¹")
# distância ao átomo mais próximo (para localizar o erro de reconstrução)
idx = np.stack(np.meshgrid(*[np.arange(k) for k in n], indexing="ij"), -1).reshape(-1, 3)
rpts = org + (idx / n) @ A; dmin = np.full(len(rpts), 1e9)
for t0 in tau:
    for sh in np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).reshape(3, -1).T:
        dmin = np.minimum(dmin, np.linalg.norm(rpts - t0 - sh @ A, axis=1))
qq = np.linspace(0, shells[keep][-1], 6001); rr = np.linspace(0, 20, 2001)
EXTS = ("const", "quad", "g0")
def deconv(v, name):
    """ṽ por casca; v_at(r) para as 3 extrapolações abaixo da 1ª casca útil:
    const = ṽ(q1); quad = par em q pelas 2 primeiras cascas; g0 = ṽ(0) = Ω⟨v⟩/nat (média absoluta do cubo), par em q."""
    Vg = np.fft.fftn(v) / N
    vt = np.array([(Om * np.sum(np.conj(Sg[sel]) * Vg[sel]) / np.sum(np.abs(Sg[sel]) ** 2)).real
                   for sel in (ins & (key == g) for g in shells)])
    qk, vk = shells[keep], vt[keep]
    def vq(x, ext):
        r = np.interp(x, qk, vk, right=0.0); lo = x < qk[0]
        if ext == "quad": r[lo] = vk[0] + (vk[1] - vk[0]) / (qk[1] ** 2 - qk[0] ** 2) * (x[lo] ** 2 - qk[0] ** 2)
        elif ext == "g0": v0q = Om * v.mean() / nat; r[lo] = v0q + (vk[0] - v0q) * x[lo] ** 2 / qk[0] ** 2
        else: r[lo] = vk[0]
        return r
    j0 = np.sinc(np.outer(rr, qq) / np.pi)
    VR = {e: np.trapezoid(vq(qq, e)[None, :] * qq ** 2 * j0, qq, axis=1) / (2 * np.pi ** 2) for e in EXTS}
    # reconstrução Σ_{G≠0} ṽ S/Ω e^{iGr} contra v limitado à esfera (G = 0 fora)
    vtG = np.zeros(n); vtG[ins] = vq(Gn[ins], "const")
    err = (np.fft.ifftn(np.where(ins, Vg, 0) * N).real - np.fft.ifftn(vtG * Sg / Om * N).real).reshape(-1)
    vb = np.fft.ifftn(np.where(ins, Vg, 0) * N).real.reshape(-1)
    say(f"# [{name}] ṽ(q1={qk[0]:.2f}) = {vk[0]:+.2f}, Ω⟨v⟩/nat = {Om * v.mean() / nat:+.2f} Ha·bohr³; reconstrução: rms(v) {vb.std() * HA:.3f},"
        f" erro rms {np.sqrt(np.mean(err ** 2)) * HA:.3f}, máx {np.abs(err).max() * HA:.3f} eV (em d_átomo = {dmin[np.argmax(np.abs(err))]:.2f} bohr)")
    say(f"#   erro rms por d_átomo (bohr) 0–1/1–2/2–3/3–4/4+: " + "/".join(
        f"{np.sqrt(np.mean(err[(dmin >= a) & (dmin < b)] ** 2)) * HA:.3f}" for a, b in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 99))))
    say(f"#   v_at(r) eV, const|quad|g0: " + " ".join(f"r={r}:" + "|".join(f"{VR[e][r * 100] * HA:+.2f}" for e in EXTS) for r in (2, 3, 4, 5, 6, 7, 8, 10)))
    return VR, vt
VRp = {}; vtab = {"q": shells, "keep": keep, "r": rr}
for k, p in parts.items():
    VRp[k], vtab[f"vt_{k}"] = deconv(p, k)
    for e in EXTS: vtab[f"vat_{k}_{e}"] = VRp[k][e]
np.savez(D + f"a11_vat_{PFX}.npz", **vtab)

# ---------- 2. PAO e elementos de matriz ----------
P = open(PAO).read().split("\n")
def block(tag):
    i0 = [k for k, l in enumerate(P) if l.strip().startswith("<" + tag)][0]
    i1 = [k for k, l in enumerate(P) if l.strip().startswith(tag + ">")][0]
    return np.array([l.split() for l in P[i0 + 1:i1]], float)
R0 = block("pseudo.atomic.orbitals.L=0"); R1 = block("pseudo.atomic.orbitals.L=1")
rp = R0[:, 1]; Rs = R0[:, 2:5].T; Rpz = R1[:, 2:5].T; rcut = rp[-1]
def orbitals(X, Y, Z):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2); rs = np.where(r > 1e-12, r, 1e-12)
    phis = [np.interp(r, rp, Rs[z], right=0) / np.sqrt(4 * np.pi) for z in range(3)]
    c = np.sqrt(3 / (4 * np.pi))
    for z in range(3):
        R = np.interp(r, rp, Rpz[z], right=0)
        phis += [c * R * X / rs, c * R * Y / rs, c * R * Z / rs]   # px, py, pz (ordem OpenMX, conferida em diag_sp)
    return np.array(phis)                                          # s1 s2 s3 p1x p1y p1z p2x ... (= scfout)
def matel(dk, h, vrs, chunk=40):
    """S e −Σ_k⟨μ|v(|r−d_k|)|ν⟩ (eV) para cada potencial radial de vrs (dict nome → v(rr))."""
    x = np.arange(-rcut, rcut + h / 2, h)
    S = np.zeros((12, 12)); V = {k: np.zeros((12, 12)) for k in vrs}
    for i0 in range(0, len(x), chunk):
        X, Y, Z = (a.ravel() for a in np.meshgrid(x[i0:i0 + chunk], x, x, indexing="ij"))
        s = X ** 2 + Y ** 2 + Z ** 2 < rcut ** 2; X, Y, Z = X[s], Y[s], Z[s]
        ph = orbitals(X, Y, Z); S += ph @ ph.T * h ** 3
        dist = [np.sqrt((X - d[0]) ** 2 + (Y - d[1]) ** 2 + (Z - d[2]) ** 2) for d in dk]
        for k, vr in vrs.items():
            pot = -sum(np.interp(dd, rr, vr) for dd in dist)
            V[k] += (ph * pot) @ ph.T * h ** 3 * HA
    return S, V

dat = np.load(D + "a11_blocos.npz")
SP = list(range(12)); s1p1 = [0, 3, 4, 5]
def metrics(dH, S):
    c = np.sum(S[np.ix_(s1p1, s1p1)] * dH[np.ix_(s1p1, s1p1)]) / np.sum(S[np.ix_(s1p1, s1p1)] ** 2)
    pb = dH[3:6, 3:6]; ap = pb - np.trace(pb) / 3 * np.eye(3); w, U = np.linalg.eigh(ap)
    k = np.argmax(np.abs(w)); oth = [j for j in range(3) if j != k]
    return c, w[k] - w[oth].mean(), U[:, k], dH[0, 3:6], np.linalg.norm(dH[np.ix_([0, 1, 2], list(range(3, 12)))])
H = 0.10   # convergido: h = 0,25…0,10 bohr iguais em 4 dígitos (rodada 1)
allv = {f"{p}|{e}": VRp[p][e] for p in parts for e in EXTS}
say(f"# h = {H} bohr, raio r_c = {rcut:.2f} bohr. escalar = tr[SΔH]/tr[S²] em {{s1,p1}}; σ−π = λ_eixo − média(outros), parte sem traço do bloco p1;"
    " v = ΔH[s1,(px,py,pz)1]; ‖sp‖ = Frobenius s(3ζ)–p(3ζ); cos(eixo), cos(v) contra a DFT. eV.")
termo = {}
for i in (0, 1, 2, 30, 31, 32):
    dk = dat[f"mis{i}"] / B2A
    dHd = (dat[f"Hs{i}"] - dat[f"Hb{i}"])[np.ix_(SP, SP)]; Sd = dat[f"S{i}"][np.ix_(SP, SP)]
    md = metrics(dHd, Sd)
    if len(dk):
        Sn, V = matel(dk, H, allv)
    else:
        Sn, V = Sd, {k: np.zeros((12, 12)) for k in allv}
    if i == 1: say(f"# |S_PAO − S_scfout|máx (átomo 1) = {np.abs(Sn - Sd).max():.1e}")
    for e in EXTS: V[f"total|{e}"] = sum(V[f"{p}|{e}"] for p in parts)
    termo[i] = (V, dHd, Sd)
    say(f"# át {i}, {len(dk)} ausentes | fonte          | escalar  σ−π    eixo p1            | v_s1p1 x y z              | |v|    | ‖sp‖   | cos(eixo) cos(v)")
    rows = [("DFT", md)] + [(k, metrics(V[k], Sd)) for k in ([f"total|{e}" for e in EXTS] + [f"{p}|{e}" for e in EXTS for p in parts])]
    for k, (c, sp_, axv, vec, nsp) in rows:
        if len(dk) == 0 and k != "DFT" and k != "total|const": continue
        ca = abs(axv @ md[2]); cv = vec @ md[3] / (np.linalg.norm(vec) * np.linalg.norm(md[3]) + 1e-30)
        say(f"{i:3d} | {k:14s} | {c:+.4f} {sp_:+.4f} {np.round(axv, 2)!s:18s} | {vec[0]:+.4f} {vec[1]:+.4f} {vec[2]:+.4f} | {np.linalg.norm(vec):.4f} | {nsp:.4f} | {ca:.2f} {cv:+.2f}")
for e in EXTS:
    say(f"# resto ΔH^DFT − ΔH^aus (total, ext {e}): át | escalar σ−π | |v| | ‖sp‖")
    for i in (0, 1, 30, 32):
        V, dHd, Sd = termo[i]; c, sp_, axv, vec, nsp = metrics(dHd - V[f"total|{e}"], Sd)
        say(f"{i:3d} | {c:+.4f} {sp_:+.4f} | {np.linalg.norm(vec):.4f} | {nsp:.4f}")
np.savez(D + f"a11_termo_{PFX}.npz", **{f"dH_{i}_{k.replace('|', '_')}": V[k] for i, (V, _, _) in termo.items() for k in V},
         **{f"dHdft_{i}": termo[i][1] for i in termo})
out.close()
