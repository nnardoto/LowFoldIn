# Surface spectral analysis — contracts (v1)

Authoritative interface specification for `atomselect`, `spectral`, `specmap`
and their modules. Implementation must follow this file; any deviation must be
reported and decided by the orchestrator, who then updates this file.

## §0. Conventions

- All indices are 1-based. Atom order = the order of the structure file, which
  must equal the HSIO atom order. Atoms are never reordered.
- Units: energies in **eV**; lengths in **Å** (Cartesian); k-points in
  **fractional reciprocal coordinates** of the cell described by the HSIO;
  accumulated k-distance in **1/Å** computed with the reciprocal metric
  `b_i · a_j = 2π δ_ij`.
- Kinds: `dp` from module `types`; on-disk integers `int32`/`int64` and reals
  `real64` from `iso_fortran_env`.
- Errors in library routines: return `ierr /= 0` plus `errmsg` (character(len=*)
  or allocatable) whenever stated; executables print the message and `stop 1`.
  Library routines must not `stop` unless stated.
- New source files carry the SPDX/copyright header used by the project.
- Fourier convention (existing, `core::compute_Hk_at_kpoint`):
  `H(k) = Σ_R H(R) exp(+i 2π k·n_R)`, lattice gauge (no τ phases). Projections
  onto atoms are gauge-invariant (diagonal phase per orbital), so this is valid.

## §1. Geometry input (decision for C0.1)

HSIO stores no geometry, species or lattice. **The HSIO format is not changed.**
Geometry comes from the existing structure file format already read by
`unfold::ReadStructureFile`:

```
&STRUCTURE
    LatticeVectors(1,:) = a1x, a1y, a1z    ! rows = a1,a2,a3, Å
    LatticeVectors(2,:) = ...
    LatticeVectors(3,:) = ...
/
&SPECIES            ! optional
<nSpecies>
<label> <nOrb>      ! nOrb counted in the HSIO basis (includes spinor doubling)
/
&ATOMS
<nAtoms>
<label> x y z       ! Cartesian Å, same order as the HSIO
/
```

Lines starting with `!` or `#` and blank lines are skipped inside the
`&SPECIES`/`&ATOMS` blocks. The per-atom `label` is the species label; it is
used as given (no element symbols are invented). If `&SPECIES` is present the
`spectral` program checks `nOrb(label(a)) == HSIO%Orbs(a)` for every atom.

### Module `structure_io` (src/mod_structure_io.f90)

```fortran
integer, parameter, public :: LABEL_LEN = 8
type, public :: structure_t
    integer :: nat = 0
    real(dp) :: lattice(3,3) = 0     ! lattice(i,:) = a_i, Å
    real(dp), allocatable :: cart(:,:)          ! (3,nat) Å
    character(len=LABEL_LEN), allocatable :: label(:)   ! (nat)
    logical :: has_species = .false.
    integer :: nspecies = 0
    character(len=LABEL_LEN), allocatable :: species_name(:)
    integer, allocatable :: species_norb(:)
end type
subroutine read_structure(filename, st, ierr, errmsg)       ! never stops
function   structure_fingerprint(st) result(fp)            ! character(len=16)
function   fnv1a64_hex(text) result(hex)                    ! character(len=16)
function   reciprocal_lattice(st) result(b)                 ! real(dp) (3,3), b(i,:) = b_i, 1/Å, incl. 2π
function   cart_to_frac(st, r) result(f)                    ! real(dp) (3)
function   species_norb_of(st, iat) result(n)               ! -1 if no species table
```

`unfold::ReadStructureFile` keeps its public signature and behaviour; it may be
refactored into a thin wrapper around `read_structure` (error → same message
prefix and `stop 1`).

### Fingerprint (§1.1)

Canonical text, lines terminated by `achar(10)`:
```
LFSTRUCT1
<nat>
<label_1> <ix_1> <iy_1> <iz_1>
...
```
`label` = `trim(adjustl(label))`; `ix = nint(x * 1e4_dp, int64)` written with
`I0` (so values are rounded to 1e-4 Å and `-0` never appears). Lattice vectors
are NOT included (atoms identity only).
Hash: **FNV-1a 64-bit** over the bytes of the canonical text
(offset basis 0xcbf29ce484222325, prime 0x100000001b3, arithmetic mod 2^64),
output as 16 lowercase hex chars. Implementation must avoid signed-integer
overflow (e.g. 32-bit limbs held in int64). Test vectors:
`fnv1a64("") = cbf29ce484222325`, `fnv1a64("a") = af63dc4c8601ec8c`,
`fnv1a64("foobar") = 85944171f73967e8`.

## §2. `regions.nml`

```fortran
&REGION_FILE
    format_version = 1
    natoms         = 12
    fingerprint    = '0123456789abcdef'     ! optional ('' = none)
    species        = 'Ga', 'N', ...         ! optional; nat labels in order
/
&REGION
    name  = 'top'
    atoms = 1, 2, 3, 4, 5, 6
/
&REGION
    name  = 'bottom'
    atoms = 7, 8
/
```

- `&REGION_FILE` must be first and appear once. `format_version` must be 1.
- Any number of `&REGION` blocks; read sequentially until `iostat_end`.
- Buffers: `atoms` and `species` are read into buffers of size
  `MAX_REGION_ATOMS = 100000` initialised to a sentinel (`-huge(0)` / `''`).
  The atom count is the number of leading non-sentinel entries; a non-sentinel
  value after a sentinel is an error ("gap in atoms list").
- Region `name`: 1–64 chars from `[A-Za-z0-9_.]`, must not be `total`
  (reserved), unique (case-sensitive).
- Validation errors (all fatal, message names region and offending value):
  index < 1, index > natoms, duplicated index inside a region, empty region,
  duplicated name, invalid name, natoms ≠ structure nat, species list ≠
  structure labels, fingerprint ≠ structure fingerprint, namelist syntax error
  (report iomsg).
- Fingerprint absent (`''`): accepted only when the caller passes
  `allow_no_fingerprint = .true.`, and then a warning is printed.
- Regions may overlap. Stored atom lists keep the user order (sorted
  ascending by atomselect when it writes them).
- Partition: `is_partition = (every atom in exactly one region)`.
- Species labels longer than `LABEL_LEN` (8) are truncated by the namelist
  read (limitation; labels written by `atomselect` always fit).

### Module `atomic_regions` (src/mod_atomic_regions.f90)

```fortran
integer, parameter, public :: REGION_NAME_LEN = 64, MAX_REGION_ATOMS = 100000
type, public :: atom_region_t
    character(len=REGION_NAME_LEN) :: name = ''
    integer, allocatable :: atoms(:)
end type
type, public :: region_set_t
    integer :: format_version = 1
    integer :: natoms = 0
    character(len=16) :: fingerprint = ''
    logical :: has_species = .false.
    character(len=LABEL_LEN), allocatable :: species(:)
    integer :: nregion = 0
    type(atom_region_t), allocatable :: regions(:)
end type
subroutine read_regions_file(filename, rs, ierr, errmsg)         ! syntax + intrinsic checks
subroutine write_regions_file(filename, rs, ierr, errmsg)        ! overwrites; canonical layout
subroutine validate_regions(rs, st, allow_no_fingerprint, ierr, errmsg)  ! vs structure
subroutine add_region(rs, name, atoms, ierr, errmsg)             ! name/duplicate/range checks
logical function regions_are_partition(rs)
logical function is_valid_region_name(name)
```

### Module `atomic_selection` (src/mod_atomic_selection.f90)

```fortran
subroutine parse_index_list(text, nmax, list, ierr, errmsg)
! text: "1-12", "15, 17, 22", "1-10, 15-20" (spaces optional, commas separate
! items, ranges inclusive a-b with a<=b). list: allocatable, sorted ascending.
! Errors: empty text, non-integer token, reversed range, index <1 or >nmax,
! duplicates (reported, not silently removed).
```
Geometric selection and layer routines (§6) are added to this module in F8.

## §3. `spectral.nml`

```fortran
&SPECTRAL
    h_file          = 'slab_H.hsio'      ! spin_mode 'none' or 'noncollinear'
    h_up_file       = ''                 ! spin_mode 'collinear'
    h_dn_file       = ''                 ! spin_mode 'collinear'
    s_file          = 'slab_S.hsio'      ! required
    structure_file  = 'slab.nml'         ! required (§1)
    regions_file    = 'regions.nml'      ! optional; '' → only 'total'
    allow_no_fingerprint = .false.
    spin_mode       = 'none'             ! 'none' | 'collinear' | 'noncollinear'
    method          = 'mulliken'         ! 'mulliken' | 'lowdin'
    energy_reference = 'unknown'         ! 'fermi' | 'absolute' | 'unknown'
    fermi_energy_ev = 0.0                ! metadata only, used iff fermi_known
    fermi_known     = .false.
    output_prefix   = 'slab'             ! writes <prefix>.spectral.bin
    write_text      = .false.            ! also <prefix>.spectral.txt
    tol_norm          = 1.0e-8
    tol_negative_warn = 1.0e-3
    tol_s_min         = 1.0e-10
/
&KPATH
    nvertex = 3
    vertex(:,1) = 0.3333333333, 0.3333333333, 0.0
    vertex(:,2) = 0.0, 0.0, 0.0
    vertex(:,3) = 0.5, 0.0, 0.0
    label   = 'K', 'G', 'M'
    npoints = 40, 30                     ! intervals per segment (nvertex-1 values)
/
```

- **Energy reference.** `spectral` NEVER shifts energies. Eigenvalues are
  written exactly as produced by the HSIO Hamiltonian (eV).
  `openmx2hsio` shifts H by −E_F·S by default, so for default conversions the
  user should set `energy_reference='fermi'` (then E_F = 0 by definition,
  `fermi_known` is forced true with value 0). For `--no-fermi-shift` files use
  `'absolute'` and optionally `fermi_known=.true., fermi_energy_ev=<E_F in eV>`.
  `'unknown'` is the safe default; maps are then labelled "E (eV, reference
  unknown)".
- **Spin.**
  - `none`: one H. nspin = 1.
  - `collinear`: `h_up_file`, `h_dn_file` share `s_file`; nspin = 2 independent
    channels (channel 1 = up, 2 = down).
  - `noncollinear`: one spinor H; nspin = 1. HSIO `Orbs(a)` already counts both
    spin components and the real ordering is atom-major
    `[TNO_a up | TNO_a down]` (`ReorderSpinorAtomMajor`). Region→orbital
    mapping uses the HSIO offsets `off(a) = Σ_{b<a} Orbs(b)` and takes all
    `Orbs(a)` orbitals of atom a, so all spin components are included without
    assuming any internal spin ordering.
- H and S HSIO files must agree in nAtoms, NOrbs, Orbs, Neighbors, IndexList;
  nAtoms must equal the structure nat.
- k-path (`kpath` module, §3.1): `MAX_KVERTEX = 100`.

### §3.1 Module `kpath` (src/mod_kpath.f90)

```fortran
integer, parameter, public :: KLABEL_LEN = 16
type, public :: kpath_t
    integer :: nk = 0, nvertex = 0
    real(dp), allocatable :: kfrac(:,:)     ! (3,nk)
    real(dp), allocatable :: kdist(:)       ! (nk) 1/Å, kdist(1)=0
    integer,  allocatable :: vertex_k(:)    ! (nvertex) index of vertex in k list
    character(len=KLABEL_LEN), allocatable :: label(:)   ! (nvertex)
end type
subroutine build_kpath(vertex, npoints, label, recip, path, ierr, errmsg)
! vertex(3,nvertex) fractional; npoints(nvertex-1) >= 1; recip(3,3) rows b_i.
! Segment s: points vertex_s + j/npoints(s)*(vertex_{s+1}-vertex_s),
! j = 0..npoints(s)-1; then the final vertex is appended once.
! nk = sum(npoints) + 1. No duplicated junction points, uniform steps.
! kdist increments = norm2( matmul(dk_frac, recip) ) (Cartesian 1/Å).
! A zero-length segment (equal consecutive vertices) is an error.
```
Not reused from `unfold::ReadPCPath` on purpose: that generator uses the
fractional norm (not the reciprocal metric), a non-uniform last step and no
labels; changing it would alter `unfold` output. Documented as a known
duplication.

### §3.2 Module `spectral_engine` (src/mod_spectral_engine.f90)

```fortran
integer, parameter, public :: MAX_KVERTEX = 100, SPEC_PATH_LEN = 512
type, public :: spectral_input_t      ! &SPECTRAL fields (§3) + &KPATH (vertex, label, npoints, nvertex)
type, public :: spectral_timing_t     ! read_s, build_s, diag_s, proj_s, write_s, total_s (measured)
subroutine read_spectral_input(filename, inp, ierr, errmsg)   ! &SPECTRAL and &KPATH, any order
subroutine run_spectral(inp, sd, timing, ierr, errmsg)
subroutine process_kpoint(n, ngroup, Hr, Sr, Rpt, kfrac, groups, Hk, Sk, C, SC, E, W, &
                          t_build, t_diag, t_proj, ierr)        ! per-k body, caller-owned workspaces
```
Energy-reference labels: `fermi` → "E - E_F (eV), E_F = 0 by definition";
`absolute` → "E (eV)" (+ E_F if known); `unknown` → "E (eV, reference unknown)".
`sd%source` = "H=<h> S=<s>" (collinear: "H_up=.. H_dn=.. S=..").
Limitation: `HSDataFormat%ReadFile` stops on corrupt HSIO files (existing
behaviour, not changed); `run_spectral` only pre-checks existence.

## §4. Projection module `spectral_projection` (src/mod_spectral_projection.f90)

```fortran
type, public :: orbital_groups_t
    integer :: ngroup = 0, norb = 0
    integer, allocatable :: start(:)   ! (ngroup+1); group g = orb(start(g):start(g+1)-1)
    integer, allocatable :: orb(:)
end type
pure subroutine atom_orbital_offsets(orbs_per_atom, offset)   ! offset(nat+1), offset(1)=0
subroutine build_orbital_groups(orbs_per_atom, nlist, list_start, list_atoms, groups)
! Generic CSR atom lists → orbital lists (group g atoms = list_atoms(list_start(g):list_start(g+1)-1)).
subroutine overlap_times_vectors(n, nband, S, C, SC)           ! SC = S·C via ZHEMM('L','U'); S untouched
subroutine mulliken_weights(n, nband, C, SC, groups, W)        ! W(nband, ngroup)
subroutine lowdin_weights(n, nband, S, C, groups, W, tol_s_min, smin, ierr)   ! F9
subroutine generalized_eigensolve(n, H, S, E, C, ierr)         ! wraps zhegv(1,'V','U') on COPIES
```

- Mulliken: `W(b,g) = Re Σ_{μ∈g} conj(C(μ,b)) · SC(μ,b)`. Negative values are
  kept as-is (no clipping, no abs). Mulliken weights are populations, not
  probabilities.
- Löwdin: `S^{1/2}` from `zheev` of a COPY of S; if min eigenvalue ≤
  `tol_s_min` → `ierr=1`, `smin` returned (caller aborts naming the k-point).
  Never floor eigenvalues. `W(b,g) = Σ_{μ∈g} |(S^{1/2}C)(μ,b)|²`.
- `generalized_eigensolve`: copies H and S, `zhegv(1,'V','U')`, eigenvalues
  ascending, `C^† S C = I`. Uses workspace query (lwork=-1). Never `stop`s;
  `ierr = info`. (Equivalent to the private
  `unfold::SolveGeneralizedEigenproblem`, which `stop`s on error and therefore
  cannot be used by a library caller that must report the k-point; unfold is
  left unchanged.)
- The S passed to projection must be the original S(k), never the zhegv output.
- `total` is a group containing all orbitals; for Mulliken it equals 1 per band
  by construction — a numerical check, not a physical result.

## §5. Spectral data file (module `spectral_data`, src/mod_spectral_data.f90)

```fortran
integer, parameter, public :: SPIN_NONE = 0, SPIN_COLLINEAR = 1, SPIN_NONCOLLINEAR = 2
integer, parameter, public :: METHOD_MULLIKEN = 1, METHOD_LOWDIN = 2
integer, parameter, public :: EREF_UNKNOWN = 0, EREF_FERMI = 1, EREF_ABSOLUTE = 2
type, public :: spectral_data_t
    integer :: version_major = 1, version_minor = 0
    integer :: nk = 0, nband = 0, nspin = 1, nregion = 0, nvertex = 0
    integer :: spin_mode = SPIN_NONE, method = METHOD_MULLIKEN, energy_ref = EREF_UNKNOWN
    logical :: fermi_known = .false.
    real(dp) :: fermi_energy_ev = 0
    character(len=16)  :: fingerprint = ''
    character(len=512) :: source = ''          ! free text: input file names
    real(dp), allocatable :: kfrac(:,:)        ! (3,nk)
    real(dp), allocatable :: kdist(:)          ! (nk) 1/Å
    integer,  allocatable :: vertex_k(:)       ! (nvertex)
    character(len=16), allocatable :: vertex_label(:)    ! (nvertex)
    character(len=64), allocatable :: region_name(:)     ! (nregion); region 1 = 'total'
    integer,  allocatable :: region_natoms(:)  ! (nregion)
    real(dp), allocatable :: energy(:,:,:)     ! (nband, nk, nspin) eV
    real(dp), allocatable :: weight(:,:,:,:)   ! (nband, nk, nregion, nspin)
end type
subroutine write_spectral_file(filename, sd, ierr, errmsg)
subroutine read_spectral_file(filename, sd, ierr, errmsg)
subroutine write_spectral_text(filename, sd, ierr, errmsg)
integer function find_region(sd, name)     ! 0 if absent
```

Binary layout (`access='stream', form='unformatted'`, native byte order with
canary; all fields contiguous, no record markers):

| # | field | type |
|---|---|---|
| 1 | magic `"LFSPECT"//achar(0)` | char(8) |
| 2 | version_major (=1), version_minor (=0) | int32 ×2 |
| 3 | endian canary `168496141` (0x0A0B0C0D) | int32 |
| 4 | sizeof int32, int64, real64 (=4,8,8) | int32 ×3 |
| 5 | nk, nband, nspin, nregion, nvertex | int32 ×5 |
| 6 | spin_mode, method, energy_ref, fermi_known(0/1) | int32 ×4 |
| 7 | fermi_energy_ev | real64 |
| 8 | fingerprint | char(16) |
| 9 | source | char(512) |
| 10 | kfrac(3,nk), kdist(nk) | real64 |
| 11 | vertex_k(nvertex) | int32 |
| 12 | vertex_label(nvertex) | char(16) each |
| 13 | region_name(nregion) | char(64) each |
| 14 | region_natoms(nregion) | int32 |
| 15 | energy(nband,nk,nspin) (Fortran order) | real64 |
| 16 | weight(nband,nk,nregion,nspin) (Fortran order) | real64 |
| 17 | n_ext (=0 in v1.0) | int32 |
| 18 | n_ext × { tag char(16), nbytes int64, payload } | — |
| 19 | trailer `"LFSPEND"//achar(0)` | char(8) |

Reader checks: magic, `version_major == 1` (any minor accepted; extension
blocks skipped by nbytes), canary, type sizes, dimensions > 0 (nvertex ≥ 0,
nregion ≥ 1, nspin ∈ {1,2}, nspin=2 iff spin_mode=collinear), region 1 named
`total`, vertex_k within 1..nk, file size ≥ bytes required before reading each
array, trailer present at the expected position (truncation/corruption →
`ierr/=0` with message). Units are fixed by the format: eV and 1/Å.

Text export `<prefix>.spectral.txt`: `#` header with all metadata, then one row
per (spin, k, band): `ispin ik kdist k1 k2 k3 iband energy w(1..nregion)`, with a
`#` line naming the region columns. For debugging only; not read back.

## §6. `atomselect` CLI (tools/atomselect.f90)

```
atomselect inspect  <structure.nml> [--xyz out.xyz] [--table out.txt]
atomselect fingerprint <structure.nml>
atomselect select   <structure.nml> --regions <regions.nml> --name NAME --atoms "LIST"
atomselect validate <structure.nml> <regions.nml> [--allow-no-fingerprint]
# F8:
atomselect layers   <structure.nml> --normal nx ny nz [--tol T] [--table out] [--regions file]
atomselect select   ... --normal nx ny nz --smin A --smax B     (geometric)
atomselect surface  <structure.nml> --regions <file> --normal nx ny nz --tol T --top N --bottom N [--top-name top] [--bottom-name bottom]
```
- `inspect` table: `Index Element X Y Z` (Å), original order. XYZ export uses
  the labels as given.
- `select` adds one region to `--regions` file: if the file does not exist it
  is created with a header (format_version, natoms, fingerprint, species from
  the structure); if it exists it is read, validated against the structure,
  the region is added (name must be unique), and the file is rewritten
  (canonical layout, previous regions preserved in order).
- Exit code 0 on success, 1 on any error (message on stdout prefixed
  `[atomselect error]`).

## §7. `specmap.nml` and specmap (tools/specmap.f90)

```fortran
&SPECMAP
    spectral_file = 'slab.spectral.bin'
    spin_channel  = 1                  ! 1..nspin
    kernel        = 'lorentzian'       ! 'lorentzian' (width = HWHM η) | 'gaussian' (width = σ)
    width_ev      = 0.05
    emin = -3.0, emax = 3.0, ne = 601  ! E_i = emin + (i-1)(emax-emin)/(ne-1)
    maps          = 'total', 'top', 'top-bottom', 'top+bottom'   ! up to 32
    occupation    = .false.            ! multiply by Fermi–Dirac
    temperature_k = 0.0
    mu_ev         = 0.0
    mu_set        = .false.            ! required true unless energy_ref = fermi (then μ=0 default)
    output_prefix = 'slab'
/
```
- Map expression: region names joined by `+`/`-` (e.g. `top-bottom`,
  `a+b-c`); linear combination with coefficients ±1 of the per-region
  intensities, computed after broadening without any normalisation; sign
  preserved. Overlapping regions in a sum may double-count atoms (documented).
- `A_P(k,E) = Σ_n w_n^P K(E-ε_n)`; Lorentzian `K = (1/π) η/(x²+η²)`, Gaussian
  `K = exp(-x²/2σ²)/(σ√(2π))`.
- Fermi–Dirac `f = 1/(1+exp((E-μ)/k_B T))`, `k_B = 8.617333262e-5 eV/K`,
  overflow-safe; T = 0 → step with f(μ) = 1/2. Applied as `I = f·A`.
- Output `<prefix>.<mapname>.map.dat`, mapname with `+`→`_plus_`,
  `-`→`_minus_`. Columns `k_distance energy_eV intensity`, one block per k
  separated by a blank line (pm3d). `#` header: source file, map expression,
  method, kernel+width, energy window, occupation (T, μ), units, energy
  reference, E_F if known.
- Also writes `<prefix>.kticks.gp` containing a gnuplot `set xtics (...)`
  command with vertex labels at their kdist, `set xrange`, and, if the Fermi
  level is known in the file's reference, `EF = <value>` (else `EF = NaN`).
- specmap uses only `spectral_data`, `spectral_broadening`, `spectral_map`
  (no HSIO, no LAPACK, no structure).

### Modules
```fortran
! spectral_broadening (src/mod_spectral_broadening.f90)
integer, parameter, public :: KERNEL_LORENTZIAN = 1, KERNEL_GAUSSIAN = 2
elemental function lorentzian(x, eta) result(y)
elemental function gaussian(x, sigma) result(y)
elemental function fermi_dirac(e, mu, temperature_k) result(f)
subroutine energy_grid(emin, emax, ne, grid)
subroutine broaden(energy, weight, grid, kernel, width, A)
! energy(nband,nk), weight(nband,nk), grid(ne) → A(ne,nk)

! spectral_map (src/mod_spectral_map.f90)
subroutine parse_map_expression(expr, sd, region_idx, coeff, ierr, errmsg)
subroutine write_map_file(filename, sd, A, grid, header_lines, ierr, errmsg)
subroutine write_kticks_file(filename, sd, ierr, errmsg)
function map_file_tag(expr) result(tag)
```

## §8. Default tolerances

`tol_norm = 1e-8`, `tol_negative_warn = 1e-3`, `tol_s_min = 1e-10`
(configurable in `spectral.nml`). Normalisation (total and, if partition,
Σ regions) violations above `tol_norm` print a warning with spin, k and band
(at most 20 lines, then a count).

## §9. Tests

fpm `[[test]]` entries with `source-dir = "test"`; shared harness
`test/mod_testing.f90` (`check`, `check_close`, `test_summary` — stops 1 on
failure). Test fixtures are generated in code or tiny files under
`test/data/`. Temporary output goes to `build/` or the system temp dir,
never committed.
