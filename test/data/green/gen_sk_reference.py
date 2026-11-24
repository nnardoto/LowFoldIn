#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Reference two-center (Slater-Koster) blocks computed WITHOUT the closed-form
Slater-Koster table.

The bond-frame matrix B (bond along +z) is diagonal in |m| (sigma, pi, delta).
For a bond direction n, a rotation R with R z = n gives

    E(n) = D_la(R) B D_lb(R)^T,

where D_l(R) is the real-harmonic rotation matrix obtained numerically by
least squares on random points (f_m(R^-1 r) = sum_m' f_m'(r) D_m'm).
For la > lb the convention of mod_sk_table is used:
E_ab(n) = (-1)^(la+lb) E_ba(n)^T with the same (sigma, pi, delta).

Orbital order (OpenMX): s; px py pz; d3z2-r2 dx2-y2 dxy dxz dyz.
Output: sk_reference.txt, read by test/test_sk_contact.f90.
"""
import numpy as np

S3 = np.sqrt(3.0)


def harmonics(l, r):
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    if l == 0:
        return np.ones_like(x)[..., None]
    if l == 1:
        return np.stack([x, y, z], axis=-1)
    if l == 2:
        return np.stack([0.5 * (3 * z * z - (x * x + y * y + z * z)),
                         0.5 * S3 * (x * x - y * y),
                         S3 * x * y, S3 * x * z, S3 * y * z], axis=-1)
    raise ValueError


def rot_matrix(l, R, rng):
    pts = rng.normal(size=(200, 3))
    f = harmonics(l, pts)                      # f_m(r)
    g = harmonics(l, pts @ R)                  # f_m(R^-1 r) (rows r^T R = (R^T r)^T)
    D, *_ = np.linalg.lstsq(f, g, rcond=None)  # g = f D
    return D


# |m| of each real harmonic in the bond frame (z = bond axis)
MABS = {0: [0], 1: [1, 1, 0], 2: [0, 2, 2, 1, 1]}
# pairing inside the bond frame: same harmonic type couples
TYPE = {0: ['s'], 1: ['x', 'y', 'z'], 2: ['z2', 'x2y2', 'xy', 'xz', 'yz']}
CANON = {'s': 'c0', 'z': 'c0', 'z2': 'c0', 'x': 'c1', 'xz': 'c1', 'y': 's1', 'yz': 's1',
         'x2y2': 'c2', 'xy': 's2'}


def bond_frame(la, lb, v):
    B = np.zeros((2 * la + 1, 2 * lb + 1))
    for a, ta in enumerate(TYPE[la]):
        for b, tb in enumerate(TYPE[lb]):
            if CANON[ta] == CANON[tb]:
                B[a, b] = v[MABS[la][a]]
    return B


def rotation_to(n, rng):
    n = n / np.linalg.norm(n)
    # random rotation about the bond axis is allowed (result must not depend on it)
    t = rng.normal(size=3)
    t -= t.dot(n) * n
    t /= np.linalg.norm(t)
    u = np.cross(n, t)
    return np.column_stack([t, u, n])          # R e_z = n


def sk_ref(la, lb, n, v, rng):
    if la > lb:
        return (-1) ** (la + lb) * sk_ref(lb, la, n, v, rng).T
    R = rotation_to(np.asarray(n, float), rng)
    return rot_matrix(la, R, rng) @ bond_frame(la, lb, v) @ rot_matrix(lb, R, rng).T


def main():
    rng = np.random.default_rng(20261119)
    cases = []
    for la in range(3):
        for lb in range(3):
            for _ in range(3):
                n = rng.normal(size=3)
                n /= np.linalg.norm(n)
                v = rng.normal(size=3)
                cases.append((la, lb, n, v))
    # special directions
    for la in range(3):
        for lb in range(3):
            for n in ([0, 0, 1], [0, 0, -1], [1, 0, 0], [0, 1, 0], [1, 1, 1]):
                n = np.array(n, float) / np.linalg.norm(n)
                cases.append((la, lb, n, np.array([0.7, -0.4, 0.2])))
    with open('sk_reference.txt', 'w') as f:
        f.write(f'{len(cases)}\n')
        for la, lb, n, v in cases:
            E = sk_ref(la, lb, n, v, rng)
            f.write(f'{la} {lb} ' + ' '.join(f'{x:.17e}' for x in n) + ' '
                    + ' '.join(f'{x:.17e}' for x in v) + '\n')
            for a in range(2 * la + 1):
                f.write(' '.join(f'{E[a, b]:.17e}' for b in range(2 * lb + 1)) + '\n')
    # sanity: p-p analytic form
    n = np.array([0.3, -0.5, 0.81]); n /= np.linalg.norm(n); v = np.array([1.3, -0.6, 0.0])
    E = sk_ref(1, 1, n, v, rng)
    A = v[0] * np.outer(n, n) + v[1] * (np.eye(3) - np.outer(n, n))
    assert np.allclose(E, A, atol=1e-12), 'p-p rotation check failed'
    print(f'wrote sk_reference.txt ({len(cases)} blocks)')


if __name__ == '__main__':
    main()
