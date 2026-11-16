# `spectral` — band energies and region projection weights along a k-path or on a k-grid

`spectral` is the **Compute** step of the surface spectral pipeline
(`atomselect` → `spectral` → `specmap`). For every k-point of a path
(`&KPATH`, for (k, E) maps) or of a 2D grid (`&KGRID`, for constant-energy
kx × ky maps) it
builds H(k) and S(k) from HSIO files, solves the generalised eigenproblem
`H(k) c = ε S(k) c`, and projects every eigenstate on the atomic regions of a
`regions.nml` file. It writes the discrete data (energies + weights) to a
binary file `<prefix>.spectral.bin`; broadening into maps is done afterwards,
without re-diagonalising, by [`specmap`](specmap.md).

```text
OpenMX .scfout ─ openmx2hsio ─> *_H.hsio, *_S.hsio ┐
structure.nml ─ atomselect ───> regions.nml        ├─> spectral ─> <prefix>.spectral.bin (+ .txt)
spectral.nml (&SPECTRAL + &KPATH|&KGRID) ──────────┘
```

## Usage

```bash
spectral <spectral.nml>
fpm run spectral -- spectral.nml
fpm run spectral --profile release -- spectral.nml      # OpenMP build
```

Exactly one argument. On error the message is printed on standard error
prefixed `[spectral error]` and the exit code is 1.

## Input: `spectral.nml`

The file contains `&SPECTRAL` and exactly one of `&KPATH` (line cut) or
`&KGRID` (2D grid), in any order. Giving both, or neither, is an error.

### `&SPECTRAL`

| field | type | default | meaning |
|---|---|---|---|
| `h_file` | string | `''` | Hamiltonian HSIO; required for `spin_mode = 'none'` or `'noncollinear'` |
| `h_up_file`, `h_dn_file` | string | `''` | spin-up / spin-down Hamiltonians; both required for `spin_mode = 'collinear'` |
| `s_file` | string | `''` | overlap HSIO (**required**; shared by both spin channels) |
| `structure_file` | string | `''` | structure file (**required**, same format as `unfold`; see [atomselect.md](atomselect.md#structure-file)) |
| `regions_file` | string | `''` | `regions.nml`; empty → only the `total` group is computed |
| `allow_no_fingerprint` | logical | `.false.` | accept a regions file with an empty fingerprint (warning) |
| `spin_mode` | string | `'none'` | `'none'`, `'collinear'` or `'noncollinear'` |
| `method` | string | `'mulliken'` | `'mulliken'` or `'lowdin'` |
| `energy_reference` | string | `'unknown'` | `'fermi'`, `'absolute'` or `'unknown'` (see below) |
| `fermi_energy_ev` | real (eV) | `0.0` | metadata only; meaningful only if `fermi_known` |
| `fermi_known` | logical | `.false.` | whether `fermi_energy_ev` is known |
| `output_prefix` | string | `'slab'` | writes `<prefix>.spectral.bin` |
| `write_text` | logical | `.false.` | also write `<prefix>.spectral.txt` |
| `tol_norm` | real | `1.0e-8` | tolerance of the normalisation checks (> 0) |
| `tol_negative_warn` | real | `1.0e-3` | Mulliken weights below `-tol_negative_warn` are reported (> 0) |
| `tol_s_min` | real | `1.0e-10` | Löwdin: abort if the smallest eigenvalue of S(k) is ≤ this (> 0) |

String values of `spin_mode`, `method` and `energy_reference` are
case-insensitive.

### `&KPATH`

| field | type | meaning |
|---|---|---|
| `nvertex` | integer | number of vertices, 2..100 |
| `vertex(:,i)` | 3 reals | vertex `i` in **fractional reciprocal coordinates** of the HSIO cell |
| `label` | strings (≤ 16 chars) | one label per vertex (e.g. `'K', 'G', 'M'`) |
| `npoints` | `nvertex-1` integers ≥ 1 | number of intervals of each segment |

Segment `s` contributes the points `vertex_s + j/npoints(s)·(vertex_{s+1} −
vertex_s)`, `j = 0..npoints(s)−1`; the final vertex is appended once. Hence
`nk = sum(npoints) + 1`, steps are uniform within a segment and **junction
points are not duplicated**. The accumulated distance `kdist` (starting at
0) uses the reciprocal metric (`b_i·a_j = 2π δ_ij`), so it is in **1/Å**
(Cartesian), not in fractional units. Two equal consecutive vertices (a
zero-length segment) are an error. A path cannot contain discontinuities
(e.g. `…K | M…`): every segment is joined to the next.

### `&KGRID` — 2D k-grid for constant-energy maps

| field | type | meaning |
|---|---|---|
| `origin` | 3 reals | corner of the grid, **fractional reciprocal coordinates** |
| `span1`, `span2` | 3 reals each | edge vectors of the grid (fractional); non-zero and not collinear |
| `n1`, `n2` | integers ≥ 2 | number of points along `span1` / `span2` (both ends included) |

Point `(i, j)` is `origin + (i−1)/(n1−1)·span1 + (j−1)/(n2−1)·span2`,
`i = 1..n1`, `j = 1..n2`, stored with the flat index `ik = i + (j−1)·n1`
(`i` fastest), so `nk = n1·n2` diagonalisations are done. The grid can lie in
any plane; for a slab with vacuum along z, use spans in the b1–b2 plane.
Examples:

```fortran
&KGRID                         ! full surface BZ, Γ at the centre (any 2D lattice)
    origin = -0.5, -0.5, 0.0
    span1  =  1.0,  0.0, 0.0
    span2  =  0.0,  1.0, 0.0
    n1 = 81, n2 = 81
/
&KGRID                         ! zoom around Γ: 20% of b1 and b2 in each direction
    origin = -0.2, -0.2, 0.0
    span1  =  0.4,  0.0, 0.0
    span2  =  0.0,  0.4, 0.0
    n1 = 101, n2 = 101
/
```

For a non-orthogonal surface cell (e.g. hexagonal) the region spanned by
`b1`, `b2` is a rhombus; the map is still plotted in Cartesian kx, ky, so it
shows the true shape. The cost grows as `n1·n2`: an 81 × 81 grid of the
12-orbital example takes ≈ 0.2 s on one thread, but for large DFT slabs the
grid is usually the dominant cost (OpenMP over k-points applies unchanged).

### Example

```fortran
&SPECTRAL
    h_file           = 'model_H.hsio'
    s_file           = 'model_S.hsio'
    structure_file   = 'model.nml'
    regions_file     = 'regions.nml'
    spin_mode        = 'none'
    method           = 'mulliken'
    energy_reference = 'fermi'
    output_prefix    = 'model'
    write_text       = .true.
/
&KPATH
    nvertex = 3
    vertex(:,1) = 0.5, 0.0, 0.0
    vertex(:,2) = 0.0, 0.0, 0.0
    vertex(:,3) = 0.5, 0.5, 0.0
    label   = 'X', 'G', 'M'
    npoints = 60, 85
/
```

## Consistency checks before the calculation

- The H and S files must agree in `nAtoms`, `NOrbs`, `Orbs`, neighbour cells
  and `IndexList`; `nAtoms` must equal the structure `nat`.
- If the structure has a `&SPECIES` block, `nOrb(label(a)) == Orbs(a)` is
  checked for every atom.
- The regions file is validated against the structure (natoms, species,
  fingerprint, indices, names) exactly as `atomselect validate` does.
- Only the existence of the HSIO files is pre-checked. A corrupt or
  truncated HSIO file stops the program inside the existing HSIO reader
  (`HSDataFormat%ReadFile`), not with a `[spectral error]` message.

## Energy reference

**`spectral` never shifts energies.** Eigenvalues are written exactly as
produced by the HSIO Hamiltonian, in eV. The HSIO files themselves store no
Fermi level, so the reference must be declared with `energy_reference`:

| situation | setting | axis label used downstream |
|---|---|---|
| `openmx2hsio` default conversion (H shifted by −E_F·S, so E_F = 0) | `energy_reference = 'fermi'` (forces `fermi_known = .true.`, `fermi_energy_ev = 0`) | "E − E_F (eV), E_F = 0 by definition" |
| conversion with `openmx2hsio --no-fermi-shift` | `energy_reference = 'absolute'`, optionally `fermi_known = .true., fermi_energy_ev = <E_F in eV>` | "E (eV)" (+ E_F if known) |
| not sure | `energy_reference = 'unknown'` (default, safe) | "E (eV, reference unknown)" |

`fermi_energy_ev` is metadata only: it is stored in the output file and used
by `specmap` for the `EF` line of the gnuplot ticks file, never to shift
energies.

## Spin

| `spin_mode` | Hamiltonians | nspin | notes |
|---|---|---|---|
| `none` | `h_file` | 1 | |
| `collinear` | `h_up_file`, `h_dn_file` (+ shared `s_file`) | 2 | two independent channels: 1 = up, 2 = down |
| `noncollinear` | `h_file` (spinor) | 1 | `Orbs(a)` already counts both spin components |

For non-collinear files `openmx2hsio` writes the orbitals in **atom-major**
order `[TNO_a up | TNO_a down]` for each atom. Region → orbital mapping uses
the HSIO offsets `off(a) = Σ_{b<a} Orbs(b)` and takes **all** `Orbs(a)`
orbitals of every atom in the region, so both spin components are always
included and no assumption on the internal spin ordering is needed.
Spin-resolved (up/down) projections inside a spinor calculation are not
provided.

## Equations

For each k-point, `C` holds the eigenvectors of `H(k)c = εS(k)c` normalised so
that `C†S(k)C = I` (LAPACK `zhegv`, eigenvalues ascending). Region `P` is the
set of orbitals μ of its atoms.

**Mulliken** (default):

```text
w_n^P(k) = Re Σ_{μ∈P} C*_{μn} (S C)_{μn}
```

These are *populations*, not probabilities: with a non-orthogonal basis they
can be **negative** (and larger than 1); they are kept as computed (no
clipping, no absolute value). Values below `-tol_negative_warn` are reported.

**Löwdin** (`method = 'lowdin'`):

```text
w_n^P(k) = Σ_{μ∈P} |(S^{1/2} C)_{μn}|²
```

`S^{1/2}` is built from the eigen-decomposition of a copy of S(k). If the
smallest eigenvalue of S(k) is ≤ `tol_s_min`, the program **aborts** naming
the spin channel, k-index and k-point; eigenvalues are never floored.
Löwdin weights are non-negative by construction.

**The `total` group.** Region 1 of the output is always `total`, a group
with all orbitals. For both methods `w^total = 1` for every band *by
construction* (for Mulliken this follows from `C†SC = I`), so it is a
**numerical check**, not a physical result.

**Normalisation checks** printed at the end:

- `max |total − 1|` (and how many bands exceed `tol_norm`);
- if the regions form a partition (every atom in exactly one region),
  `max |Σ_regions − 1|`; otherwise "not checked";
- the minimum weight (Mulliken: count below `-tol_negative_warn`; Löwdin:
  count below `-tol_norm`).

Individual violations are printed as `[spectral warning]` lines with spin,
k and band (at most 20 lines per kind, then a count).

The spectral function itself is formed by `specmap`:
`A_P(k,E) = Σ_n w_n^P(k) K(E − ε_n(k))` (see [specmap.md](specmap.md)).

## Output

### `<prefix>.spectral.bin`

Stream-access, unformatted, native byte order with an endianness canary; all
fields contiguous (no record markers). Units are fixed by the format: eV and
1/Å.

| # | field | type |
|---|---|---|
| 1 | magic `"LFSPECT"` + NUL | char(8) |
| 2 | version_major (=1), version_minor (=0) | int32 ×2 |
| 3 | endian canary `168496141` (0x0A0B0C0D) | int32 |
| 4 | sizeof int32, int64, real64 (= 4, 8, 8) | int32 ×3 |
| 5 | nk, nband, nspin, nregion, nvertex | int32 ×5 |
| 6 | spin_mode (0 none, 1 collinear, 2 noncollinear), method (1 Mulliken, 2 Löwdin), energy_ref (0 unknown, 1 fermi, 2 absolute), fermi_known (0/1) | int32 ×4 |
| 7 | fermi_energy_ev | real64 |
| 8 | fingerprint | char(16) |
| 9 | source (`"H=<h> S=<s>"`, collinear `"H_up=.. H_dn=.. S=.."`) | char(512) |
| 10 | kfrac(3,nk), kdist(nk) | real64 |
| 11 | vertex_k(nvertex) (index of each vertex in the k list) | int32 |
| 12 | vertex_label(nvertex) | char(16) each |
| 13 | region_name(nregion) (region 1 = `total`) | char(64) each |
| 14 | region_natoms(nregion) | int32 |
| 15 | energy(nband, nk, nspin), Fortran order | real64 |
| 16 | weight(nband, nk, nregion, nspin), Fortran order | real64 |
| 17 | n_ext (0 for k-path files; 1 for k-grid files) | int32 |
| 18 | n_ext × {tag char(16), nbytes int64, payload} | — |
| 19 | trailer `"LFSPEND"` + NUL | char(8) |

**k-grid files** (written with `&KGRID`) have `version_minor = 1`,
`nvertex = 0`, `kdist = 0` (not meaningful on a grid) and one extension block
`KGRID` (152 bytes): `n1, n2` (int32), `origin(3), span1(3), span2(3)` and the
reciprocal lattice `recip(3,3)` (rows b_i in 1/Å incl. 2π), all real64. With
`recip` stored, `specmap` computes Cartesian k without the structure file.
k-path files are unchanged (version 1.0, `n_ext = 0`).

The reader (used by `specmap`) checks magic, major version 1 (any minor;
unknown extension blocks skipped, `KGRID` validated), canary, type sizes, dimensions, that region 1 is
`total`, vertex indices, the file size before each array, and the trailer;
truncated or corrupt files give an error message.

### `<prefix>.spectral.txt` (if `write_text = .true.`)

For inspection/debugging only; never read back. A `#` header with all
metadata (version, dimensions, spin mode, method, energy reference, E_F,
fingerprint, source, units, vertices with `ik` and `kdist` — or, for a
k-grid, `n1`, `n2`, origin and spans —, regions with their atom counts), a `#` line naming the columns, then one row per
(spin, k, band):

```text
# ispin ik kdist k1 k2 k3 iband energy total top bottom interior
```

### Standard output

Number of OpenMP threads, calculation stages, periodic progress updates with
completed `(spin, k)` points, elapsed time and ETA, normalisation checks, a
summary (dimensions, spin mode, method, energy reference, fingerprint,
regions, files written) and a timing report. `read`, `k-loop`, `write`,
`total` are wall-clock times;
`build`, `diagonalization`, `projection` are accumulated over threads
(thread-seconds) and may exceed the k-loop wall time with several threads.

## Performance and parallelism

The loop over (spin, k) is parallelised with OpenMP when built with
`fpm build --profile release` (which adds `-fopenmp`); a plain `fpm build`
runs serially with identical results. Output is bitwise independent of the
thread count. Keep BLAS/LAPACK single-threaded inside OpenMP threads (e.g.
`OPENBLAS_NUM_THREADS=1` or `VECLIB_MAXIMUM_THREADS=1`). Eigenvectors are not
stored for all k-points; only energies and weights are kept. Details and the
measured timings (tiny 18-orbital test case only) are in
[performance.md](performance.md).

## Approximations and limitations

- Weights are projections of Kohn–Sham-like eigenstates of the given H; no
  photoemission matrix elements, photon energy/polarisation, final states,
  k_z broadening or many-body effects are included.
- k-points are either a path (line cuts) or one planar 2D grid; path
  discontinuities are not supported, and symmetry is not used to reduce the
  grid (the full `n1·n2` points are diagonalised).
- Mulliken weights can be negative and do not sum to 1 per region unless the
  regions form a partition.
- Overlapping regions are allowed; each region is projected independently,
  so summing overlapping regions double-counts atoms.
- Species labels are limited to 8 characters.
- A corrupt HSIO file stops the program in the existing HSIO reader.
- Fourier convention: `H(k) = Σ_R H(R) exp(+i2πk·n_R)` (lattice gauge, no
  intra-cell phases); atom projections are gauge-invariant, so this is valid.
- The tools have been verified on synthetic models and small HSIO fixtures
  only (see Tests); no validation against a real DFT slab is claimed.

## Band structure, weights, spectral function, maps and ARPES

- **Band structure** ε_n(k): eigenvalues of `H(k)c = εS(k)c` along the path
  (column `energy`).
- **Projection weight** w_n^P(k): how much of state n lies on region P
  (Mulliken or Löwdin; Mulliken values are populations and may be
  negative).
- **Projected spectral function** A_P(k,E) = Σ_n w_n^P(k) K(E − ε_n(k)): the
  weights broadened by an *artificial* kernel K chosen by the user (computed
  by `specmap`). It is not a many-body spectral function.
- **Intensity map**: A_P, or a ± combination of regions, optionally times a
  Fermi–Dirac occupation, on a (k, E) grid — the files written by `specmap`.
- **Experimental ARPES**: a measured photoemission intensity, which also
  depends on matrix elements, photon energy and polarisation, final states,
  k_z, resolution and many-body effects. None of these are modelled; maps
  can guide the interpretation of ARPES but are not simulations of it.

## Tests

- `fpm test test_kpath` — k-path construction (point count, no duplicated
  junctions, uniform steps, reciprocal-metric distances, errors) and the 2D
  k-grid (corners, index order, Cartesian coordinates, in-plane basis,
  errors).
- `fpm test test_projection` — Mulliken and Löwdin weights, `S·C`, the
  generalised eigensolver, the `tol_s_min` abort (synthetic matrices).
- `fpm test test_spectral_data` — binary format round trip and reader
  checks (truncation/corruption), text export.
- `fpm test test_spectral_engine` — input parsing, consistency on the
  diatomic fixture, Löwdin, partition sums, non-collinear and collinear
  spin, error cases, file round trip and bitwise determinism (1 vs 4
  threads in an OpenMP build); `&KGRID` parsing (exclusive with `&KPATH`) and
  k-grid energies against an independent `zhegv` at every grid point.
- `fpm test test_structure_io`, `test_regions` — structure and regions input.

## Worked example

[`examples/surface_spectral`](../examples/surface_spectral/README.md) runs the
full pipeline on a synthetic 8-layer non-orthogonal tight-binding slab (not
DFT data), including the `spectral.nml` shown above, with measured
normalisation checks and the energies/weights of its two surface states.
