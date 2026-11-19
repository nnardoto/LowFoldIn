#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
"""Gerador dos dados de referência do módulo de Green (fase 0 do plano v4).

Protótipo em Python/NumPy, independente do código Fortran, para uma cadeia
não ortogonal complexa (4 orbitais por camada principal, S != I, S01 != 0):

  nonorth_chain.txt   blocos H00, H01, S00, S01; para cada energia, a Green de
                      superfície g_s nos dois lados (decimação de López-Sancho
                      na forma K = zS - H), conferida contra uma cadeia finita
                      longa resolvida por recursão de Dyson;
  aa_reference.txt    DOS de bulk por PL (Mulliken com termos cruzados) por
                      integração em k_s e por embedding de dois lados (teste
                      A/A), e o valor incompleto só com Tr(G_mm S_mm).

Uso (requer numpy):  python3 gen_reference.py   (grava no diretório do script)
Os arquivos gerados são versionados; o CI não executa este script.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261119)
n = 4
eta = 0.05
energies = [-1.5, 0.2, 2.5]


def crand(shape):
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


def herm(a):
    return 0.5 * (a + a.conj().T)


H00 = herm(crand((n, n))) * 0.6 + np.diag(np.linspace(-1.0, 1.0, n))
H01 = crand((n, n)) * 0.35
S00 = np.eye(n) + herm(crand((n, n))) * 0.04
S01 = crand((n, n)) * 0.03
H10, S10 = H01.conj().T, S01.conj().T


def kblocks(e):
    z = e + 1j * eta
    return z * S00 - H00, z * S01 - H01, z * S10 - H10


def decimation(k00, k01, k10, side, tol=1e-14, maxit=200):
    # side 'top': bulk nas PLs <= 0 -> U = K10, L = K01
    u, l = (k10, k01) if side == 'top' else (k01, k10)
    ds, db = k00.copy(), k00.copy()
    for it in range(1, maxit + 1):
        gb_l = np.linalg.solve(db, l)
        gb_u = np.linalg.solve(db, u)
        ugl = u @ gb_l
        ds = ds - ugl
        db = db - ugl - l @ gb_u
        u, l = -u @ gb_u, -l @ gb_l
        if max(abs(u).max(), abs(l).max()) < tol:
            break
    return np.linalg.inv(ds), it


def finite_chain(k00, k01, k10, side, nlayer=4000):
    u, l = (k10, k01) if side == 'top' else (k01, k10)
    g = np.linalg.inv(k00)
    for _ in range(nlayer):
        g = np.linalg.inv(k00 - u @ g @ l)
    return g


def write_mat(f, a):
    for j in range(a.shape[1]):
        for i in range(a.shape[0]):
            f.write(f"{a[i, j].real: .17e} {a[i, j].imag: .17e}\n")


maxdiff = 0.0
results = []
for e in energies:
    k00, k01, k10 = kblocks(e)
    row = []
    for side in ('top', 'bottom'):
        gs, it = decimation(k00, k01, k10, side)
        gf = finite_chain(k00, k01, k10, side)
        maxdiff = max(maxdiff, abs(gs - gf).max())
        row.append(gs)
    results.append(row)

with open(os.path.join(HERE, 'nonorth_chain.txt'), 'w') as f:
    f.write('# Non-orthogonal complex chain: reference surface Green functions (gen_reference.py)\n')
    f.write(f'# max |g_decimation - g_finite_chain(4000 layers)| = {maxdiff:.3e}\n')
    f.write(f'{n} {len(energies)} {eta:.17e}\n')
    for a in (H00, H01, S00, S01):
        write_mat(f, a)
    for e, (gtop, gbot) in zip(energies, results):
        f.write(f'{e:.17e}\n')
        write_mat(f, gtop)
        write_mat(f, gbot)

# --- teste A/A: DOS de bulk por PL -------------------------------------------
nk3 = 20000
aa = []
for e in energies:
    z = e + 1j * eta
    acc = 0.0
    for k3 in (np.arange(nk3) + 0.5) / nk3:
        ph = np.exp(2j * np.pi * k3)
        hk = H00 + H01 * ph + H10 / ph
        sk = S00 + S01 * ph + S10 / ph
        acc += -np.trace(np.linalg.solve(z * sk - hk, sk)).imag / np.pi
    rho_k = acc / nk3
    k00, k01, k10 = kblocks(e)
    gA, _ = decimation(k00, k01, k10, 'top')     # meio abaixo
    gB, _ = decimation(k00, k01, k10, 'bottom')  # meio acima
    sig = k10 @ gA @ k01 + k01 @ gB @ k10
    g00 = np.linalg.inv(k00 - sig)
    g0m = -g00 @ k10 @ gA        # G_{0,-1}
    g0p = -g00 @ k01 @ gB        # G_{0,+1}
    rho_emb = -np.trace(g00 @ S00 + g0m @ S01 + g0p @ S10).imag / np.pi
    rho_diag = -np.trace(g00 @ S00).imag / np.pi
    aa.append((e, rho_k, rho_emb, rho_diag))

with open(os.path.join(HERE, 'aa_reference.txt'), 'w') as f:
    f.write('# Bulk DOS per principal layer of nonorth_chain (A/A test), gen_reference.py\n')
    f.write(f'# k_s integration with {nk3} points; columns: E rho_kint rho_embedding rho_diag_only\n')
    f.write(f'{len(aa)}\n')
    for e, a, b, c in aa:
        f.write(f'{e:.17e} {a:.17e} {b:.17e} {c:.17e}\n')

print('max |decimation - finite chain| =', maxdiff)
for row in aa:
    print('E=%5.2f  k-int %.12f  embed %.12f  diag-only %.12f  (diag error %.2f%%)'
          % (row[0], row[1], row[2], row[3], 100 * abs(row[3] - row[2]) / abs(row[2])))
