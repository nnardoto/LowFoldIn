#!/usr/bin/env python3
"""A25: termos de face por átomo central de um filme de TeBulk_fo2 (hsiocell, L8) sobre o volume semi-infinito, só do volume:
 (a) VNA+xc dos vizinhos ausentes (gerador do A12, a12_termo.py: termo(), RMAX 14,5 bohr; ausentes = sítios do volume com
     n2 >= 8 + células removidas); átomos cujo conjunto de ausentes não contém célula removida recebem o bloco do A12 do filme x1
     (a12_blocos.npz, mesma subrede e camada; idêntico por translação);
 (b) L=2 excedente do campo cristalino do volume (A20 'cfvol'): −Σ_aus β_camada (d̂d̂ − I/3) (2 camadas, ajuste do a20_cf_vol.py,
     ausentes < 4 Å) menos a parte sem traço do bloco p1 de (a); só entra se máx|R| > 5 meV (no x1 isso dá os átomos 22 e 24 do A20);
 (c) barreira FD por distância assinada: s(r) = d_presentes(r) − d_ausentes(r), V_b = V0/(1 + exp(−(s − s0)/ws)), V0 = 2 eV;
     (s0, ws) calibrados uma vez na face plana x1 contra a barreira planar do A18/A20 (V0 2, d 0,5, w 0,35), Frobenius nos
     blocos on-site 29×29 dos 12 átomos centrais; malha cartesiana h = 0,15 bohr na esfera de 7 bohr (como a18_barreira.py).
Raiz dos dados: $A25_ROOT (padrão: ~/calc/tese-Te/green_surface/), com a mesma árvore lfi_tese/... (espelho no scicom).
uso: a25_face.py calib                → calibra (s0, ws); grava trincheira_vol/calib.json
     a25_face.py check                → refaz o termo (a) de um átomo do x1 com termo() e compara com a12_blocos.npz
     a25_face.py x1|tr out.sitepot    → sitepot completo (a)+(b)+(c) da face plana x1 ou da trincheira (x11, sem a célula 5,7,0)"""
import os, sys, json, time, numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import least_squares
G = os.environ.get("A25_ROOT", "/home/nathanael/calc/tese-Te/green_surface/").rstrip("/") + "/"
D = G + "lfi_tese/"; T = D + "campo_semiinf/termo/"; TV = D + "trincheira_vol/"
B2A = 0.529177210903; RC = 7.0; H = 0.15
V0, D_PL, W_PL = 2.0, 0.5, 0.35
RMAX_A = (7.0 + 7.0 + 0.5) * B2A          # raio do A12 (14,5 bohr) em Å
L8 = 8

def readstruct(f):
    L = open(f).read().split("\n")
    lat = np.array([[float(x) for x in l.split("=")[1].split()] for l in L if "LatticeVectors" in l])
    i = [k for k, l in enumerate(L) if l.strip() == "&ATOMS"][0]; nat = int(L[i + 1])
    return lat, np.array([[float(x) for x in L[i + 2 + a].split()[1:4]] for a in range(nat)])
latb, posb = readstruct(T + "hs/TeBulk.structure.nml")
FILMS = {"x1": (T + "hs/TeFilm_x1_L8.structure.nml", 1, set()),
         "tr": (TV + "hs/TeFilm_x11_L8_tr.structure.nml", 11, {(5, 7)})}

class Film:
    def __init__(s, tag):
        f, s.n1, s.rem = FILMS[tag]; s.tag = tag
        s.lat, s.pos = readstruct(f)
        s.lab = []
        for p in s.pos:
            for b in range(3):
                fr = np.linalg.solve(latb.T, p - posb[b])
                if np.abs(fr - np.round(fr)).max() < 1e-6:
                    t = np.round(fr).astype(int); s.lab.append((b, t[0], t[1], t[2])); break
            else: raise RuntimeError(f"átomo fora da rede do volume: {p}")
        s.lab = np.array(s.lab)
        assert (s.lab[:, 3] == 0).all() and s.lab[:, 2].min() == 0 and s.lab[:, 2].max() == L8 - 1
        s.cent = np.where(s.lab[:, 2] >= 4)[0]
    def present(s, t1, t2):
        rem = np.zeros(np.shape(t1), bool)
        for (r1, r2) in s.rem: rem |= (np.mod(t1, s.n1) == r1) & (t2 == r2)
        return (t2 <= L8 - 1) & ~rem
    def sites(s, i, rad):
        """sítios do volume a < rad (Å) do átomo i: posições, máscara 'presente'."""
        b, t1, t2, _ = s.lab[i]
        n = int(np.ceil(rad / 4.0)) + 2
        g = np.stack(np.meshgrid(np.arange(t1 - n - 3, t1 + n + 4), np.arange(t2 - n, t2 + n + 1), np.arange(-n, n + 1), np.arange(3),
                                 indexing="ij"), -1).reshape(-1, 4)
        P = posb[g[:, 3]] + g[:, :3] @ latb
        m = np.linalg.norm(P - s.pos[i], axis=1) < rad
        g, P = g[m], P[m]
        return P, s.present(g[:, 0], g[:, 1]), g

# ---------- (c) barreira ----------
import subprocess, tempfile
def orbitals(rel):   # pao_eval do LowFoldIn (paochk/paoeval_bin), como a18_barreira.py
    with tempfile.TemporaryDirectory() as td:
        fi, fo = os.path.join(td, "in.bin"), os.path.join(td, "out.bin")
        with open(fi, "wb") as f:
            np.array([len(rel)], np.int32).tofile(f); np.ascontiguousarray(rel, np.float64).tofile(f)
        subprocess.run([T + "paochk/paoeval_bin", fi, fo], check=True)
        return np.fromfile(fo, np.float64).reshape(len(rel), 29).T
_GR = {}
def grid():
    if "g" not in _GR:
        x = np.arange(-RC, RC + H / 2, H); X, Y, Z = np.meshgrid(x, x, x, indexing="ij"); m = X**2 + Y**2 + Z**2 <= RC**2
        P = np.stack([X[m], Y[m], Z[m]], 1)
        _GR["g"] = (P, orbitals(P), H**3)
    return _GR["g"]
def sfield(F, i):
    """s(r) = d_pres − d_aus (Å) nos pontos da malha do átomo i; None se não há ausente a < 7 bohr + 8 Å."""
    P, phi, h3 = grid()
    S, pr, _ = F.sites(i, RC * B2A + 8.0)
    if pr.all(): return None
    r = F.pos[i] + P * B2A
    if LAM <= 0:
        dp, _ = cKDTree(S[pr]).query(r); da, _ = cKDTree(S[~pr]).query(r)
        return dp - da
    return softmin(r, S[pr]) - softmin(r, S[~pr])
LAM = 0.0     # λ (Å) da distância suave; 0 = distância ao sítio mais próximo (calib.json decide)
def softmin(r, Q, ch=20000):
    """distância suave −λ ln Σ_j exp(−|r − Q_j|/λ) (→ min_j |r − Q_j| quando λ → 0)."""
    out = np.empty(len(r))
    for a in range(0, len(r), ch):
        d = np.linalg.norm(r[a:a + ch, None, :] - Q[None], axis=2); m = d.min(1)
        out[a:a + ch] = m - LAM * np.log(np.exp(-(d - m[:, None]) / LAM).sum(1))
    return out
def bar_block(s, s0, ws):
    P, phi, h3 = grid()
    v = V0 / (1 + np.exp(-np.clip((s - s0) / ws, -700, 700)))
    M = (phi * v) @ phi.T * h3; return 0.5 * (M + M.T)
def planar_blocks(F):
    P, phi, h3 = grid(); ytop = F.pos[F.cent, 1].max(); out = {}
    for i in F.cent:
        y = F.pos[i, 1] + P[:, 1] * B2A
        v = V0 / (1 + np.exp(-np.clip((y - (ytop + D_PL)) / W_PL, -700, 700)))
        M = (phi * v) @ phi.T * h3; out[i] = 0.5 * (M + M.T)
    return out

# ---------- (a) VNA+xc ----------
_A12 = {}
def a12():
    if not _A12:
        src = open(T + "a12_termo.py").read(); src = src[:src.index("if sys.argv[1] == \"check\"")]
        ns = {"__name__": "a12"}; exec(compile(src, "a12_termo.py", "exec"), ns); _A12.update(ns)
    return _A12
def vna_block(F, i, flat):
    S, pr, g = F.sites(i, RMAX_A)
    mis = ~pr
    if not mis.any(): return None, 0, False
    b, t1, t2, _ = F.lab[i]
    remd = mis & (g[:, 1] <= L8 - 1)            # ausentes por remoção (não acima do corte)
    if not remd.any():                           # mesmo conjunto do filme plano: bloco do A12 por translação
        a = b + 3 * t2
        return (flat[f"dH{a}"] if f"dH{a}" in flat else None), int(mis.sum()), False
    A = a12(); Ri = F.pos[i] / B2A
    Sm, Vna, Vxc, neg = A["termo"](Ri, [p / B2A for p in S[mis]])
    M = Vna + Vxc; return 0.5 * (M + M.T), int(mis.sum()), True

# ---------- (b) L=2 excedente ----------
def beta2():
    d = np.load(D + "vizinhos/a11_blocos.npz"); btv, bpos = d["btv"], d["bpos"]
    def U(v): v = v / np.linalg.norm(v); return np.outer(v, v) - np.eye(3) / 3
    def tl(M): return M - np.trace(M) / 3 * np.eye(3)
    SH = [(0, 3.0), (3.0, 4.0)]
    def neigh(b):
        out = []
        for n1 in range(-2, 3):
            for n2 in range(-2, 3):
                for n3 in range(-2, 3):
                    for j in range(3):
                        R = bpos[j] + n1 * btv[0] + n2 * btv[1] + n3 * btv[2] - bpos[b]; r = np.linalg.norm(R)
                        if 0.1 < r < 4.6: out.append(R)
        return out
    def sub(i):
        r = d["pos"][i] - d["t"]
        for b in range(3):
            f = np.linalg.solve(btv.T, r - bpos[b])
            if np.abs(f - np.round(f)).max() < 1e-3: return b
    cand = {sub(i): i for i in d["at"]}
    A, y = [], []
    for b, i in sorted(cand.items()):
        Tm = tl(d[f"Hb{i}"][3:6, 3:6]); nb = neigh(b)
        Bs = [sum((U(R) for R in nb if SH[s][0] < np.linalg.norm(R) < SH[s][1]), np.zeros((3, 3))) for s in range(2)]
        for p in [(0, 0), (1, 1), (0, 1), (0, 2), (1, 2)]: A.append([Bk[p] for Bk in Bs]); y.append(Tm[p])
    beta, *_ = np.linalg.lstsq(np.array(A), np.array(y), rcond=None)
    return beta, U, tl
def l2_excess(F, i, Mvna, beta, U, tl):
    S, pr, g = F.sites(i, 4.0)
    cf = np.zeros((3, 3))
    for R in S[~pr] - F.pos[i]:
        cf -= (beta[0] if np.linalg.norm(R) < 3.0 else beta[1]) * U(R)
    ta = tl(Mvna[3:6, 3:6]) if Mvna is not None else np.zeros((3, 3))
    return cf - ta

def write_sitepot(F, blocks, out, hdr):
    with open(out, "w") as f:
        f.write("LFSITEPOT 1\n" + "".join(f"# {h}\n" for h in hdr) + f"norb 29\nnentry {len(blocks)}\n")
        for a in sorted(blocks):
            p = F.pos[a]; f.write(f"entry 0 {a + 1} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f} block\n")
            for row in blocks[a]: f.write(" ".join(f"{x:.10f}" for x in row) + "\n")

if __name__ == "__main__":
    mode = sys.argv[1]; t0 = time.time()
    if mode == "calib":
        F = Film("x1"); ref = planar_blocks(F); allres = []
        for LAM in (0.0, 0.5, 1.0, 1.5, 2.5):
          globals()["LAM"] = LAM
          sf = {i: sfield(F, i) for i in F.cent}
          act = [i for i in F.cent if sf[i] is not None]
          nref = np.sqrt(sum(np.sum(ref[i] ** 2) for i in F.cent))
          def res(p):
              return np.concatenate([(bar_block(sf[i], p[0], p[1]) - ref[i]).ravel() for i in act] +
                                    [ref[i].ravel() for i in F.cent if i not in act])
          best = None
          for s0 in (-1.0, 0.0, 0.5, 1.0):
              for ws in (0.2, 0.4, 0.8):
                  r = least_squares(res, [s0, ws], bounds=([-5, 0.02], [5, 5]))
                  if best is None or r.cost < best.cost: best = r
          s0, ws = best.x; rr = np.sqrt(2 * best.cost) / nref
          print(f"# λ = {LAM} Å: calibração (x1, 12 átomos centrais, {len(act)} com ausentes a < 7 bohr + 8 Å): s0 = {s0:.4f} Å, ws = {ws:.4f} Å;"
                f" resíduo Frobenius relativo {rr:.4f} (‖ref‖ = {nref:.4f} eV)")
          print("# át(1-based) t2 | máx|ref| | máx|ger| | máx|ger − ref| | <s1|V|s1> ref/ger | <p1y|V|p1y> ref/ger (eV)")
          for i in F.cent:
              M = bar_block(sf[i], s0, ws) if sf[i] is not None else np.zeros((29, 29))
              print(f"{i + 1:3d} {F.lab[i, 2]} | {np.abs(ref[i]).max():.4f} | {np.abs(M).max():.4f} | {np.abs(M - ref[i]).max():.4f} |"
                    f" {ref[i][0, 0]:+.4f}/{M[0, 0]:+.4f} | {ref[i][4, 4]:+.4f}/{M[4, 4]:+.4f}")
          allres.append((rr, LAM, s0, ws))
        rr, LAM, s0, ws = min(allres)
        print("# escolhido (menor resíduo):", f"λ = {LAM} Å, s0 = {s0:.4f}, ws = {ws:.4f}, resíduo {rr:.4f}")
        json.dump({"s0": s0, "ws": ws, "res_rel": rr, "V0": V0, "lam": LAM}, open(TV + "calib.json", "w"), indent=1)
        print(f"# {time.time() - t0:.0f} s"); sys.exit()
    if mode == "check":
        F = Film("x1"); flat = np.load(T + "a12_blocos.npz"); A = a12()
        for i in (21, 23):
            S, pr, g = F.sites(i, RMAX_A)
            Sm, Vna, Vxc, _ = A["termo"](F.pos[i] / B2A, [p / B2A for p in S[~pr]]); M = Vna + Vxc; M = 0.5 * (M + M.T)
            print(f"át {i + 1}: {int((~pr).sum())} ausentes; máx|M − a12_blocos| = {np.abs(M - flat[f'dH{i}']).max():.2e} eV (máx|M| {np.abs(M).max():.3f})")
        print(f"# {time.time() - t0:.0f} s"); sys.exit()
    out = sys.argv[2]; F = Film(mode); cal = json.load(open(TV + "calib.json")); LAM = cal.get("lam", 0.0)
    flat = np.load(T + "a12_blocos.npz"); beta, U, tl = beta2()
    print(f"# A25 {mode}: {len(F.pos)} átomos, {len(F.cent)} centrais; s0 = {cal['s0']:.4f} Å, ws = {cal['ws']:.4f} Å; β = {np.round(beta, 4)} eV")
    print("# át(1-based) b t1 t2 | nmis(a) novo? máx|(a)| | máx|R L2| entra | máx|(c)| | máx|total|")
    blocks = {}; nnew = 0
    for i in F.cent:
        Ma, nm, new = vna_block(F, i, flat); nnew += new
        R = l2_excess(F, i, Ma, beta, U, tl); useR = np.abs(R).max() > 5e-3
        sf = sfield(F, i); Mc = bar_block(sf, cal["s0"], cal["ws"]) if sf is not None else None
        if Mc is not None and np.abs(Mc).max() <= 1e-3: Mc = None
        M = np.zeros((29, 29))
        if Ma is not None: M += Ma
        if useR: M[3:6, 3:6] += R
        if Mc is not None: M += Mc
        if Ma is not None or useR or Mc is not None: blocks[i] = M
        b, t1, t2, _ = F.lab[i]
        if i in blocks:
            print(f"{i + 1:3d} {b} {t1:2d} {t2} | {nm:3d} {'N' if new else '-'} {np.abs(Ma).max() if Ma is not None else 0:.4f} |"
                  f" {np.abs(R).max():.4f} {'S' if useR else '-'} | {np.abs(Mc).max() if Mc is not None else 0:.4f} | {np.abs(M).max():.4f}", flush=True)
    np.savez(out.replace(".sitepot", "_blocos.npz"), **{f"M{a}": M for a, M in blocks.items()}, lab=F.lab, pos=F.pos)
    write_sitepot(F, blocks, out, [f"A25 {mode}: (a) VNA+xc A12 ({nnew} átomos recalculados) + (b) L=2 excedente cfvol (β {beta[0]:.4f} {beta[1]:.4f}, > 5 meV)"
                                   f" + (c) barreira FD s = d_pres − d_aus (λ {LAM} A), V0 {V0}, s0 {cal['s0']:.6f}, ws {cal['ws']:.6f} A; eV"])
    print(f"escrito {out}: {len(blocks)} entradas ({nnew} termos (a) novos); {time.time() - t0:.0f} s")
