# Slater–Koster contact and parameter scans (`surfgreen` model `sk`, `skfit`)

The `sk` contact model puts the A/B coupling in **by hand**, the way spin–orbit
coupling is added to a tight-binding model: a few parameters with the correct
symmetry, added to two bulks that come from DFT. It works for principal layers
(PLs) of different sizes, which the `scaled` model cannot handle, and it does not
need a DFT calculation of the junction.

What it can and cannot claim: a hand-made contact shows **mechanisms** (what the
spectrum *can* do for a given coupling, gap, registry and band alignment), not
the behaviour of a specific interface. Parameters calibrated on the bulks
(`skfit`) are a transferability hypothesis that a junction DFT can test later.

```text
OpenMX A.scfout ─ openmx2hsio --no-fermi-shift --write-structure ─> A_H/S.hsio, A.nml ─ skfit (&SKFIT) ─> A.sk ┐
OpenMX B.scfout ─ openmx2hsio --no-fermi-shift --write-structure ─> B_H/S.hsio, B.nml ─ skfit (&SKFIT) ─> B.sk ┤
                                                                                     skfit (&SKCOMBINE) ─> contact.sk
surfgreen (&INTERFACE model='sk', sk_file='contact.sk', &SCAN) ─> PREFIX_sNNN.green.bin + PREFIX.scan.txt
```

## The contact

For every atom `i` of the lower PL and `j` of the upper PL (translated by
`R_par = n1 a1 + n2 a2`) closer than `sk_cutoff_ang`, and every pair of radial
shells (`a` on `i`, `b` on `j`) that has a line in the `.sk` file:

```
V_c(d)  = lambda_h * s_c * V_c(d0) * exp(-beta (|d| - d0)),     c = sigma, pi, delta
H_AB(k_par) = Σ_{i,j,R_par} E_ab(d; V) exp(2πi (k1 n1 + k2 n2)),     S_AB = 0
```

with `d = r_j + R_par − r_i` and `E_ab` the Slater–Koster two-centre table
(s, p, d). The geometry is the real one: positions from the structure files,
the stacking vector of the lower material and `&FILM interface_offset` (extra
gap and lateral registry). The same lines are used for A below B and for B
below A (`sandwich`, `side = 'bottom'`).

- **Real and spin-diagonal**: `V ⊗ I₂` on non-collinear (spinor) atoms. The SOC of
  each material stays inside its own bulk blocks.
- **`S_AB = 0`**: orbitals on the two sides are treated as orthogonal; `S` stays
  positive definite whatever the parameters.
- **f shells** are never coupled.
- `lambda_h → 0` gives exactly the `none` model (superposition without interaction).

### Conventions

- Orbital order inside an atom (OpenMX): for `l = s, p, d, f`, for each radial
  function (zeta), the `m` components `s; px, py, pz; d3z²−r², dx²−y², dxy, dxz, dyz`.
  Non-collinear atoms are `[all spin-up | all spin-down]` (`openmx2hsio`, atom-major).
- A basis is the OpenMX string: `sk_basis = 'Te s3p2d2', 'Sb s3p2d2'`. The number
  of orbitals per atom in the `.hsio` must be `nspatial` or `2·nspatial`.
- A parameter always refers to the SK table with the **lower `l` first**; the
  species of each shell says on which atom it sits. With this convention Harrison's
  universal values have the usual signs (`ssσ < 0`, `spσ > 0`, `ppσ > 0`, `ppπ < 0`).
- The table was checked against an independent reference (bond-frame matrix rotated
  with numerically built real-harmonic rotation matrices, `test/data/green/gen_sk_reference.py`).

### `.sk` file

```text
# speciesA shellA  speciesB shellB   d0(Å)  beta(1/Å)  V_sigma  V_pi  V_delta (eV)
Te p1  Sb p1   3.10  1.20   1.35  -0.42  0.00
Te s1  Sb p1   3.10  1.20   0.80
```

One line per pair of shells (`s1` = first radial s function, `p2` = second p, …),
valid in both directions. Missing trailing channels are 0; channels beyond
`min(la, lb)` are ignored. Comments start with `#`.

## `surfgreen` input

```fortran
&INTERFACE
    model = 'sk'
    sk_file = 'contact.sk'
    sk_cutoff_ang = 6.0                        ! pair cutoff (Å)
    sk_basis = 'Bi s3p2d2', 'Te s3p2d2', 'Sb s3p2d2'
    lambda_h = 1.0                             ! global factor
    sk_scale_sigma = 1.0, sk_scale_pi = 1.0, sk_scale_delta = 1.0
    sk_beta = -1.0                             ! >= 0 replaces beta of every line
/
```

The gap and registry are set in `&FILM`: `interface_offset = dx, dy, dz` (Å, added
to the natural stacking position of the film).

## Parameter scans (`&SCAN`)

```fortran
&SCAN
    parameter = 'lambda_h'
    values = 0.0, 0.25, 0.5, 1.0, 2.0          ! explicit values ...
    ! n = 21, vmin = 0.0, vmax = 2.0           ! ... or n evenly spaced values
/
```

| `parameter` | changes | contact model |
|---|---|---|
| `lambda_h` | global factor of `H_AB` | `scaled`, `sk` |
| `lambda_s` | factor of `S_AB` | `scaled` |
| `sk_scale_sigma`, `sk_scale_pi`, `sk_scale_delta` | one channel | `sk` |
| `sk_beta` | decay of all lines (1/Å) | `sk` |
| `film_shift_ev` | band alignment `H_B + Δ S_B` (also of a B lead) | any |
| `offset_x`, `offset_y`, `offset_z` | one Cartesian component of `interface_offset` (Å) | any |
| `offset_normal` | extra distance along the substrate normal (Å), added to `interface_offset` | any |

One run writes `PREFIX_s001.green.bin … PREFIX_sNNN.green.bin` and the index
`PREFIX.scan.txt` (`index value file`). Each file carries `scan_parameter`,
`scan_value`, `scan_index` and the contact parameters in its metadata
(`surfgreen --info`). The semi-infinite substrate is decimated **once per
(k, E)** and reused for every value; only a scan that shifts a semi-infinite
medium (`film_shift_ev` in `junction` mode) recomputes that medium. The output
is bitwise identical for any number of OpenMP threads.

## `skfit`: calibrating on the bulks

### `&SKFIT`: fit to the coupling between PLs of one bulk

With one PL per stacking cell the block `H_{0,1}` is the coupling across the
stacking gap (for Bi₂Te₃ with one quintuple layer per PL, the Te–Te coupling
across the van der Waals gap). `skfit` fits SK lines to it:

- **Target, anchored at E0:** `T = H_{0,1} − E0 S_{0,1}`. Replacing `(H01, S01)` by
  `(T, 0)` is exact at `E = E0` and errs by `(E − E0) S01` away from it. That is
  the same approximation the contact makes (`S_AB = 0`).
- **Spinors:** the spatial target is `Re[(T↑↑ + T↓↓)/2]`; what is dropped is reported.
- **Model:** linear least squares (LAPACK `dgels`) in the `V_c` for a given global
  `beta`, which is found by a coarse + fine grid search, or fixed. `d0` of a line
  is the shortest distance among the pairs it describes.

```fortran
&SKFIT
    h_file = 'Bi2Te3_H.hsio', s_file = 'Bi2Te3_S.hsio', structure_file = 'Bi2Te3.nml'
    stacking_axis = 3, pl_cells = 0
    basis = 'Bi s3p2d2', 'Te s3p2d2'
    fit_shells = 'Bi s1 p1', 'Te s1 p1'      ! default: first zeta of s, p, d
    cutoff_ang = 6.0
    anchor_ev = ...                          ! default: E_F from the .hsio metadata
    beta = -1.0                              ! < 0: grid search in [beta_min, beta_max]
    beta_min = 0.0, beta_max = 4.0
    output = 'Bi2Te3.sk'
    check_kpar = 0.0, 0.0, check_nk = 41, check_window_ev = 1.0
/
```

The report gives:

- the residual on the fitted shells;
- the relative error on the **whole** `H01 − E0 S01` (all shells, all `R_par`);
- the dropped spin/imaginary part;
- a **band test**: bulk eigenvalues along `k_s` at `check_kpar`, within
  `±check_window_ev` of E0, compared with those of the anchored coupling and of
  the SK coupling (max and rms, meV).

If SK cannot reproduce the bulk's own stacking gap, it is not reliable at the
interface either.

### `&SKCOMBINE`: contact between species that were never neighbours

```fortran
&SKCOMBINE
    file_a = 'Bi2Te3.sk', file_b = 'Sb.sk'
    species_a = 'Te'                ! species of A facing the film
    species_b = 'Sb'
    output = 'contact.sk'
/
```

For each X of `species_a` and Y of `species_b`, and each pair of shells present in
both the homonuclear lines X–X (from `file_a`) and Y–Y (from `file_b`):

```
V_XY = sgn · sqrt(|V_XX · V_YY|)
```

Each value is first brought to the mean `d0` by its own exponential law;
`beta = mean`. Opposite signs give 0 and are reported. This is a starting point
for the scan, not a prediction.

## Suggested workflow for Bi₂Te₃/Sb

1. Bulks in surface-adapted cells, Sb strained to the Bi₂Te₃ in-plane lattice:
   `openmx2hsio X.scfout --no-fermi-shift --write-structure`.
2. `skfit` on each bulk. Check the band test before going on.
3. `skfit` `&SKCOMBINE` → `contact.sk`.
4. `surfgreen`, mode `film`:
   - `model = 'none'` gives the reference without interaction;
   - `model = 'sk'` with `&SCAN lambda_h` (strength), `offset_normal` (gap) and
     `film_shift_ev` (alignment);
   - compare `ΔA = A(λ) − A(0)` per layer.
5. When a junction DFT exists, compare with it (the `extracted` model, phase 4).

## Validation

| test | checks |
|---|---|
| `test_sk_table` | SK blocks (s, p, d) = independent rotation reference to 1e-12; parity rules; bond along z; basis strings |
| `test_sk_contact` | `skfit` recovers `beta` and every `V` (s, p, d, heteronuclear) of a synthetic SK material; the junction, sandwich and film with the `sk` contact equal the exact `H01` (`scaled`, λ = 1) to 1e-10 with spinors and a tilted stacking vector, coherent maps included; `lambda_h = 0` equals `none`; `&SCAN` of `lambda_h`, `film_shift_ev` and `offset_normal` equals separate runs; causality; input errors |
| `examples/sk_contact` | full CLI pipeline on synthetic materials with PLs of different sizes |

## Limitations

- `S_AB = 0`, and one global `beta` per fit.
- Spin-diagonal coupling only: an interface Rashba term is not included.
- f shells are not coupled.
- The combination rule is an assumption.
- The contact does not change the on-site blocks next to the interface (no charge
  transfer); use `film_shift_ev` and `layer_potential_ev` for electrostatics.
