# `surfgreen` — Green functions of surfaces, films and junctions

`surfgreen` computes the retarded Green function of a stack of principal
layers (PLs) coupled to semi-infinite media, in the **non-orthogonal** LCAO
basis of OpenMX (or the orthonormal basis of Wannier90), and writes
layer/region-resolved spectral functions `A_r(k_par, E)` and coherent
`q_z` maps to a versioned binary file `<prefix>.green.bin`.

It implements phases 0–3 of the LowFoldIn-Green plan (v4, Nov. 2026):
PL blocks, semi-infinite decimation, central region (RGF), projections,
contact models `none` and `scaled`. The `extracted` contact (blocks taken
from a junction DFT, phase 4) is not implemented yet.

```text
OpenMX .scfout ─ openmx2hsio --no-fermi-shift --write-structure ─> A_H/S.hsio, A.nml, B_H/S.hsio, B.nml
                                                                          │
surfgreen --write-structure in.nml ─> <prefix>_central.nml ─ atomselect ─> regions.nml
                                                                          │
surfgreen in.nml ─> <prefix>.green.bin ─ surfgreen --slice ─> *.map.dat / *.cemap.dat ─> tools/gnuplot
```

The result is a one-particle Kohn–Sham spectrum of frozen `H` and `S`.
Region maps are projected spectral functions, not ARPES simulations, and the
contact model sets what a junction can show: a `scaled` contact shows a
mechanism (what *can* happen for a given `λ`, `Δ`), not the behaviour of a
specific interface.

## Usage

```bash
surfgreen INPUT.nml                    # calculation
surfgreen --write-structure INPUT.nml  # only <prefix>_central.nml (for atomselect)
surfgreen --info FILE.green.bin        # metadata, regions, diagnostics
surfgreen --slice FILE.green.bin [--map total|REGION|coh|incoh] [--iq N] [--energy E] [--out PREFIX]
```

`--slice` exports text compatible with the existing gnuplot scripts:

| k sampling | files | script |
|---|---|---|
| `&KPATH` | `PREFIX.map.dat` (`k_dist E intensity`), `PREFIX.kticks.gp` | `tools/gnuplot/map_positive.gp` |
| `&KGRID` | `PREFIX.cemap.dat` (`kx ky kz k_u k_v intensity`) at the energy closest to `--energy`, `PREFIX.cegrid.gp` | `tools/gnuplot/cemap.gp` |

Run with OpenMP threads over `(k_par, E)` points and a single-threaded BLAS
(`OPENBLAS_NUM_THREADS=1`). The output file is bitwise identical for any
number of threads.

The BLAS/LAPACK library must be safe to call from several OpenMP threads at
once. Some OpenBLAS builds (pthreads, without locking) return wrong results
in that case — the test suite catches it (`fpm test --profile release` with
`OMP_NUM_THREADS > 1`). Linking the OpenMP variant (`libopenblaso`) or
reference LAPACK fixes it. The decimation is memory-bandwidth bound: on a
16-core/32-thread CPU, one thread per physical core
(`OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=spread`) was about 35 %
faster than 32 threads with SMT.

## Formulation

All algebra uses the resolvent `K(z) = zS − H`, `z = E + iη`, with LU
factorisations and linear solves only.

**Principal layers.** The bulk is given in a surface-adapted cell: two
in-plane lattice vectors and a stacking vector `a_s` (`stacking_axis`). A PL
groups `P = pl_cells` cells; with `R = n_s a_s + R_par` the blocks are

```
H_{m,m+p}(k_par) = Σ_{R_par} H_p(R_par) exp(2πi(k1 n1 + k2 n2)),   p ∈ {−1, 0, 1}
```

and the same for `S`. `|n_s| ≤ P` is checked: hoppings are never truncated
(error with the minimal `pl_cells`; `pl_cells = 0` chooses it).

**Semi-infinite media** (López-Sancho decimation, non-orthogonal form). With
`D = K00`, `U`, `L` the couplings towards the bulk (`L ≠ U†` for `η > 0`):

```
D_s ← D_s − U g_b L,  D_b ← D_b − U g_b L − L g_b U,  U ← −U g_b U,  L ← −L g_b L,  g_b = D_b⁻¹
```

until `max(|U|,|L|) < decim_tol`. The residual `‖(K00 − U g_s L) g_s − I‖`
and the iteration count are stored per point.

**Embedding and RGF.** `Σ_A = K_CA g_A K_AC`; the central region is
block-tridiagonal and RGF gives `G_mm`, `G_{m,m±1}` and the boundary blocks
`G_{1,A} = −G_11 K_{1,A} g_A`.

**Layer DOS in a non-orthogonal basis.** Mulliken (default) keeps the cross
terms with neighbouring layers and with the semi-infinite media:

```
ρ_m(E, k_par) = −(1/π) Im Σ_{n=m−1}^{m+1} Tr[G_mn S_nm]
```

Using only `Tr(G_mm S_mm)` breaks the A/A test (1–2 % error in the reference
model). Regions are sums of per-orbital values `d_i = −(1/π) Im Σ_j G_ij S_ji`.
`projection = 'lowdin'` uses `S_C^{1/2} G_C S_C^{1/2}` restricted to the
central region (full `G_C`, cost `O(N³)` per point).

**Energies.** A rigid shift of a material is `H → H + Δ S` (`shift_ev`),
applied also to its semi-infinite medium. A layer potential `V_m`
(`layer_potential_ev`, film layers) enters as `H_mn → H_mn + (V_m + V_n)/2 S_mn`.

**Coherent map.** For a covariant detection vector
`w_i = a_i exp(i q_z d_i) exp(−d_i / 2λ)` (`d_i` = depth below the topmost
atom, `a_i` = 1 on `&COHERENT region`, `λ = escape_length_ang`), in the
atomic gauge `U_ii = exp(i k_par · r_i)` (which corrects the layer phase of a
tilted `a_s`):

```
I_coh = −(1/π) Im w'† G w',  w' = U w      (solved from K x = w', no full G)
I_incoh = −(1/π) Im Σ_m w'_m† G_mm w'_m
```

`I_coh ≥ 0` for any `w`. Varying `q_z` changes weights, never pole positions.

## Stacking

Material A is the substrate, B the film; layers are listed bottom → top
along `+a_s`:

| `mode` | `side = 'top'` | leads (bottom / top) |
|---|---|---|
| `surface` | `[A × central_layers]` | A / vacuum |
| `film` | `[A × central_layers][B × nlayers]` | A / vacuum |
| `junction` | `[A × central_layers][B × nlayers]` | A / B |
| `sandwich` | `[A × c][B × nlayers][A × c]` | A / A |

`side = 'bottom'` mirrors the stack (free surface facing `−a_s`). Between
PLs of the same material the coupling is the material's own `H_{m,m+1}`;
between A and B it is the contact block:

| `&INTERFACE model` | `H_AB`, `S_AB` |
|---|---|
| `none` | 0 (decoupled) |
| `scaled` | `λ_H H_{m,m+1}`, `λ_S S_{m,m+1}` of the `reference` material (`'substrate'` or `'film'`); PL sizes must match |
| `extracted` | phase 4 — not implemented |

The film structure must be strained to the substrate's in-plane lattice
(checked to `lattice_tol_ang`); k-points are fractional coordinates of the
substrate's reciprocal lattice.

## Input

```fortran
&GREEN
    mode = 'film'              ! 'surface' | 'film' | 'junction' | 'sandwich'
    eta_ev = 0.002
    emin = -1.0, emax = 1.0, ne = 401
    decim_tol = 1.0e-10
    decim_maxit = 200
    projection = 'mulliken'    ! or 'lowdin'
    output_prefix = 'sb_on_bi2te3'
    regions_file = ''          ! regions.nml of <prefix>_central.nml; '' = one region per PL
    allow_no_fingerprint = .true.
    vacuum_ang = 20.0          ! vacuum added to a_s of the central structure
    lattice_tol_ang = 0.05
/
&SUBSTRATE
    h_file = 'Bi2Te3_H.hsio', s_file = 'Bi2Te3_S.hsio'   ! s_file = '' -> S = I (Wannier90)
    structure_file = 'Bi2Te3_hex.nml'
    stacking_axis = 3, pl_cells = 1, side = 'top'
    shift_ev = 0.0
    central_layers = 0         ! substrate PLs inside the central region (default 1 in 'surface', else 0)
/
&FILM
    h_file = 'Sb_H.hsio', s_file = 'Sb_S.hsio'
    structure_file = 'Sb_strained.nml'
    nlayers = 6, pl_cells = 1
    shift_ev = -0.15
    layer_potential_ev = 0.0, 0.0, ...   ! optional V_m per film layer
    interface_offset = 0.0, 0.0, 0.0     ! Å, geometry only (regions, depths)
/
&INTERFACE
    model = 'scaled'           ! 'none' | 'scaled' | 'extracted' (phase 4)
    reference = 'substrate'
    lambda_h = 1.0, lambda_s = 1.0
/
&COHERENT                      ! optional
    nqz = 21, qz_min = 0.0, qz_max = 2.0   ! 1/Å
    escape_length_ang = 6.0                ! 0 = no attenuation
    region = 'film_top'                    ! '' = all central orbitals
/
&KGRID
    origin = -0.2, -0.2, 0.0
    span1 = 0.4, 0.0, 0.0
    span2 = 0.0, 0.4, 0.0
    n1 = 81, n2 = 81
/
```

`&KPATH` (as in `spectral`) replaces `&KGRID` for `E × k` cuts. The
stacking component of every k-vector must be zero.

## Energy reference and heterostructures

`openmx2hsio` shifts `H` by `−E_F S` by default. That aligns the Fermi levels
of two bulks, which is wrong for a heterostructure. Convert both materials
with

```bash
openmx2hsio Bi2Te3.scfout --no-fermi-shift --write-structure
```

and set the band alignment explicitly through `shift_ev` (from a junction or
slab calculation, core levels, …). `openmx2hsio` now stores the absolute
`E_F`, whether the shift was applied, and the energy reference in an optional
`meta` dataset of the `.hsio` container (old readers ignore it). `surfgreen`
copies them to the output metadata (`energy_reference`,
`substrate_fermi_ev`, `film_fermi_ev`) and warns when a Fermi-shifted file is
used in a heterostructure.

## Output: `.green.bin` (v1)

Versioned container in the style of `.hsio` v2 (header, datasets, TOC; no
HDF5):

| dataset | content |
|---|---|
| `meta` | `key = value` lines: mode, η, tolerances, shifts, contact model, input fingerprints, projection and gauge conventions, energy reference |
| `axes` | k-points (fractional, Cartesian in 1/Å, path distance), energies (eV), `q_z` (1/Å), path labels |
| `maps` | `total(E,k)` and `A_r(E,k)` per region |
| `coherent` | `I_coh`, `I_incoh (E,k,q_z)` (only with `&COHERENT`) |
| `diagnostics` | decimation iterations (negative = not converged) and residual per point, `min −Im G_ii`, smallest eigenvalue of `S_C` per k |

## Validation (gates)

| gate | test | checks |
|---|---|---|
| 1 | `test_layer_blocks` | `H00 + H01 e^{2πik_s} + H10 e^{−2πik_s} = H(k)` to 1e-12 (synthetic and example `.hsio`, three stacking axes); `|n_s| > pl_cells` is an error; PL regrouping preserves bands; atomic gauge with tilted `a_s` |
| 2 | `test_surface_gf` | analytic chain to 1e-10; non-orthogonal complex chain equal to the Python reference to 1e-12; residual; `−Im g_ii ≥ 0` |
| 2 | `test/test_surfgreen_cli.sh` | output bitwise identical for 1, 2 and 4 threads; `--info`, `--slice`, `--write-structure`, error exits |
| 3 | `test_central_green` | RGF = dense inverse; A/A recovers the bulk DOS to 1e-8 (diagonal-only term does not); thick slab = semi-infinite; decoupled film = finite-slab eigenvalues; rigid shift exact; coherent poles independent of `w`; Löwdin/Mulliken trace |
| 3 | `test_green_engine` | the same properties through the full input → `.green.bin` path, including A/A in `sandwich` and `junction` |

The reference data in `test/data/green/` come from
`test/data/green/gen_reference.py` (NumPy prototype; not run by CI).

## Cost

The substrate decimation dominates: `O(n³)` per `(k_par, E)` for a PL of `n`
orbitals, with iterations growing like `log(1/η)` and near band edges (see
the `niter` diagnostic). The central region adds a few `O(n³)` operations per
layer; each extra `q_z` value is one block-tridiagonal solve.

## Example

`examples/surfgreen_model/run_example.sh` runs the full pipeline on
synthetic models (NOT DFT): substrate surface, film on substrate with regions
from `atomselect`, a constant-energy map and a coherent `q_z` scan.
