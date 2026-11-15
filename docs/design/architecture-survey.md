# LowFoldIn architecture survey (verified against code, branch feat/surface-spectral-analysis @ 692012e)

## 1. Layout and fpm.toml
- `src/`: mod_types, mod_hs_io, mod_core, mod_config, mod_bands_io, mod_terminal_ui, mod_unfold, mod_openmx_convert, mod_scfout_io. They are compiled into libLowFoldIn.a.
- `tools/`: openmx2hsio.f90, unfold.f90. `app/main.f90` (program LowFoldIn). `examples/`. No `test/`, no `docs/` before this file.
- `fpm.toml`: `auto-executables/auto-tests/auto-examples = false`; `link = ["lapack","blas"]`; implicit typing and implicit-external are off.
- Executables: `openmx2hsio` (tools/openmx2hsio.f90) and `unfold` (tools/unfold.f90) ONLY. `app/main.f90` is NOT BUILT (it has no `[[executable]]` entry).
- `[features]`: standard (-std=f2023), openmp (-fopenmp), warnings, strict-warnings, -O0..-O3, debug-symbols, runtime-checks, fpe traps, sanitizers, march variants, lto, fast-math.
- Profiles: debug, debug-sanitize, development, balanced, release, portable, native, ryzen7950x, ryzen7950x-fast. All of them include `openmp`. A plain `fpm build` without `--profile` gets NO -fopenmp, so the `!$` lines are inactive.
- `fpm build` (default profile): **compiled successfully** (both executables).

## 2. Shared types
- `src/mod_types.f90` module `types`: kinds `sp=real32`, `dp=real64`, `qp=real128`. No other types.
- `src/mod_hs_io.f90` module `HS_io`, `type HSDataFormat` (l.98-113):
  - `Title character(len=256)`, `Neighbors integer` (= number of R cells), `nAtoms integer`, `NOrbs integer`
  - `Degenerecency(:) integer, allocatable` [Neighbors] (always 1 from the converter; NOT used in the Fourier sum)
  - `Orbs(:) integer, allocatable` [nAtoms] = orbitals per atom
  - `IndexList(:,:) integer, allocatable` [3, Neighbors] = integer cell indices (n1,n2,n3)
  - `Data(:,:,:) complex(dp), allocatable` [NOrbs, NOrbs, Neighbors], always dense in memory
  - TBPs: `ReadFile`, `WriteFile`, `DisplayInfo`, `Clear`. A private `TOCEntry` type.
- Each file holds ONE matrix (H or S). H and S are separate files.

## 3. HSIO v2 binary layout (mod_hs_io.f90)
- Unformatted stream, native endianness.
- Header: `HSIO_MAGIC="LFIOFMT"//char(0)` (8 B), header_bytes i32, ver_major=2 i32, ver_minor=0 i32, endian canary 0x0A0B0C0D i32, toc_offset i64 (1-based POS), toc_count i32 (always 1).
- Payload (`write_blocksparse_hs_payload` l.414, schema `SCHEMA_BLOCKSPARSE_HS_V1=1`): Title, Neighbors, nAtoms, NOrbs, Degenerecency(Neighbors), Orbs(nAtoms), IndexList(3,Neighbors), nBlocks i64, BlockIndex(3,nBlocks) i32 (a1,a2,Rn), BlockOffset(nBlocks+1) i64, BlockData complex128. Inside a block the first index varies fastest. All-zero atom-pair blocks are dropped.
- TOC entry: name char(64) ("H"/"S", picked by the `is_overlap_like` heuristic), kind (`DATASET_BLOCKSPARSE_HS=1`), schema_ver, offset i64, nbytes i64. The reader takes the first entry with kind=1 and skips the others. The multi-dataset container exists but is unused.
- `ReadFile` rebuilds the dense `Data` with zeros filled in. The reader stops on bad magic, wrong major version, or wrong endianness.
- Metadata present / absent:
  - number of atoms: `nAtoms`. Orbitals per atom: `Orbs(:)`. Total: `NOrbs`.
  - orbital ordering: atom-contiguous blocks in atom order. Offsets are recomputed by each consumer. There is no per-orbital l/m label. **NOT PRESENT.**
  - species / element symbols: **NOT PRESENT**.
  - atomic positions (any units): **NOT PRESENT**.
  - lattice vectors: **NOT PRESENT** (only integer `IndexList`).
  - spin type flag: **NOT PRESENT**. It is inferred only from the free text in `Title` and from Orbs = 2×TNO.
  - energy units: **NOT PRESENT** as a field. By convention H is in eV from openmx2hsio (stated only in the Title text).
  - Fermi energy: **NOT PRESENT**. It is applied implicitly, because H is shifted by default.

## 4. openmx2hsio (tools/openmx2hsio.f90, src/mod_openmx_convert.f90, src/mod_scfout_io.f90)
- CLI: `openmx2hsio file.scfout [--out-prefix P] [--fermi E_F_Hartree] [--no-fermi-shift]`.
- `ReadScfoutFile` (mod_scfout_io) reads SCFOUT v3.
  - It **reads and discards** `atv` (Cartesian cell vectors), `tv`/`rtv` (lattice vectors) and `Gxyz` (positions, l.131-200). Geometry is available in the .scfout but is thrown away.
  - It keeps TotalNumOrbs, FNAN, natn, ncn, atv_ijk, Hks, iHks, OLP, ChemP, E_Temp, Solver.
- Units: `HARTREE_TO_EV=27.211386245988` is applied to H only (`Hdata = Hdata*HARTREE_TO_EV`). S is not scaled. There is no Bohr->Å conversion because no lengths are written.
- Fermi shift: `ApplyFermiShift(H,S,E_F)`: `H(R) <- H(R) - E_F*S(R)` for every R (l.422). It is ON BY DEFAULT.
  - E_F comes from ChemP in the .scfout, or from the `--fermi` override (in Hartree).
  - `--no-fermi-shift` disables it.
  - It is called for all spin cases (openmx2hsio.f90 l.154, 160, 166, 183). E_F is not stored in the file.
- Files written:
  - SpinP_switch 0: `P_S.hsio`, `P_H.hsio`.
  - SpinP_switch 1 (collinear): `P_S.hsio` (shared), `P_H_up.hsio`, `P_H_dn.hsio`, each spatial NOrbs.
  - SpinP_switch 3 (noncollinear/SOC): `P_H.hsio`, `P_S.hsio` with NOrbs = 2·Nspatial and `Orbs = 2*TotalNumOrbs`.
- Noncollinear assembly:
  - `AssembleNonCollinearH` builds a spin-major matrix [[Huu,Hud],[Hdu,Hdd]]. Hdu comes from the Hermitian mirror of Hud(-R).
  - `DoubleBlockDiagonal` builds S = diag(S,S).
  - `ReorderSpinorAtomMajor` (l.94-136) then permutes to the final order: **for each atom a, [all TNO_a orbitals spin-up, then all TNO_a orbitals spin-down]**, atom after atom. This is NOT interleaved per orbital.
  - Global index of (atom a, orbital o, spin s): `newoff(a) + o + (s==dn)*TNO_a`, where `newoff(a+1) = newoff(a) + 2*TNO_a`.
- Data convention: `Data(i,j,Rn)` = <i in cell 0 | H | j in cell IndexList(:,Rn)>. This follows from `AssembleDense`, with row = ct_AN in the home cell and col = neighbour Gh_AN in cell ncn. The cell list order is atv_ijk, and index 1 is (0,0,0).
- `HermiticityReport` prints a diagnostic of H(Γ).

## 5. unfold (tools/unfold.f90 + src/mod_unfold.f90)
- Input: one file is given as argv(1).
  - `ReadUnfoldNamelist` reads `&UNFOLD`: ProjectName, HDataFile, SDataFile, OutputDir, M(3,3), Broadening, Emin, Emax, nE, WriteSpectralFunction, PCStructureFile, SCStructureFile, MatchTolerance. It calls `ClassifyM` (diagonal check, nCells=|det M|).
  - Then the atom map is built either automatically (`BuildAtomMapFromPositions` when both structure files are given) or by hand (`ReadAtomMap`, a `&ATOMMAP` block with lines "PCAtom n1 n2 n3").
  - Then `ReadPCPath` reads the `&PCPATH` block.
- `NextDataLine` (l.173) skips blank lines and lines starting with `!`/`#`.
- `ReadStructureFile` (l.554), custom .nml format:
  - `&STRUCTURE LatticeVectors(i,:)` namelist. Rows are a1,a2,a3, in Å.
  - Optional `&SPECIES` block: `nSp`, then "El nOrb" lines.
  - `&ATOMS` block: `nAt`, then "El x y z" lines, Cartesian Å.
  - It returns FRACTIONAL positions via `(Lᵀ)⁻¹ r`, plus element symbols and the species→nOrb table. It does not return lattice vectors. Geometry exists only here, never in HSIO.
- `RunUnfolding` (l.787):
  - Builds orbit tables: AtomOrbOffset, OrbitID, OrbitMemberOrb.
  - Per PC k-point (OpenMP):
    - ScK = mod(fk·M, 1) (row vector).
    - `compute_Hk_at_kpoint` for H and S.
    - `SolveGeneralizedEigenproblem`.
    - `MatrixSqrtHermitian`, then C' = S^½ C (dense `matmul`).
    - Weights W: closed-form phases for diagonal M; projection operator with `BuildGroupTable`/`FindCellIndex` for general M.
    - SumCheck = ΣW.
- Outputs:
  - `SaveScatter` writes `<base>.unfold.dat` with columns (x, E, W).
  - `SaveSpectralFunction` writes `<base>.unfold_spectral.dat`.

## 6. H(k)/S(k) and Fourier convention
- `core::compute_Hk_at_kpoint(Hr, Rpt, Hk, kpt)` (src/mod_core.f90 l.234, `pure`, PUBLIC): `H(k) = Σ_R Data(:,:,R) · exp(+i·2π·k·n_R)`.
  - k is fractional (reduced reciprocal coordinates of the cell the HSIO describes).
  - n_R = integer IndexList.
  - **No atomic positions (τ)**: this is the lattice/periodic gauge.
  - It is dense: every R is summed, including zero blocks.
- It is the single reusable routine. It is used by `calc_band_structure`, `calc_band_structure_non_orth` and `RunUnfolding` (applied to S as well).

## 7. Solvers
- `unfold::SolveGeneralizedEigenproblem(H,S,n,Eval,Cvec)` (mod_unfold.f90 l.1256):
  - Copies H and S into Awork/Bwork, then calls `zhegv(1,'V','U',…)` (not zhegvd).
  - lwork = 2n-1 (the minimum). Stops if info ≠ 0.
  - Eigenvalues come out ascending. Eigenvectors are S-orthonormal (LAPACK itype=1: CᴴSC = I).
- `unfold::MatrixSqrtHermitian(S,n,Ssqrt,U,w)` (l.1295): copies S, calls `zheev('V','U')`, floors eigenvalues at 1e-10, returns S^½ = U·diag(√w)·Uᴴ.
- **Both routines are PRIVATE** (they are missing from the `public` list, l.59-62). Reuse needs an export or a move to a shared linear-algebra module.
- `core` declares interfaces for zheev/zheevd/zhegv/zhegvd. It uses only zheev('N') and zhegv(1,'N','U') (eigenvalues only). zhegvd and zheevd are declared but NOT USED.

## 8. Orbital → atom mapping
- There is no stored offset array. Each consumer recomputes a prefix sum of `Orbs(:)`:
  - `AtomOffset` in `DisplayInfo` and in the payload read/write routines
  - `AtomOrbOffset` in RunUnfolding (l.844)
  - `atom_offset` via `OrbitalOffsets` in mod_openmx_convert (private)
- There is no shared helper.

## 9. Spin
- Unpolarized: one H file.
- Collinear: two H files plus one S file. unfold is run once per channel. There is no spin field and no combined handling.
- Noncollinear/SOC: a spinor 2N basis in atom-major [up block, dn block] per atom, with `Orbs` doubled. unfold treats spinor orbitals as ordinary orbitals (the PC orbit identity includes spin), see examples/UnfoldSpinorNonCollinear.
- There is no spin projection, Sz, or spin-texture code: **NOT PRESENT**.

## 10. Spectral function / k-path
- `SaveSpectralFunction` (l.1390): Lorentzian kernel `A(k,E) = Σ_N W_N · (γ/π) / ((E-E_N)² + γ²)` with `γ = 0.5*Broadening`. **Broadening is the FWHM; γ is the HWHM.**
  - The energy grid is uniform over [Emin, Emax], nE points.
  - It is computed serially at write time from W, not stored.
  - Output is text triples (x, E, A) with a blank line between k-points.
- There is no Gaussian kernel and no DOS, PDOS, or orbital/atom/layer projection: **NOT PRESENT**.
- `config::BandsConfig` (mod_config.f90 l.17): npaths, `total_kpoints` (integer, allocatable scalar), nkpts(:), kpoints(3,npaths), FullPath(3,:), xpoints(:).
- `unfold::PCPathConfig` (mod_unfold.f90 l.85): the same fields plus a non-allocatable total_kpoints. It is a duplicated type.
- Path generation (`ReadBandsConfig` l.118 and `ReadPCPath` l.691, the same algorithm):
  - Each vertex line is "k1 k2 k3 nk", in fractional coordinates.
  - Segment i produces nkpts(i) points, starting at vertex i and excluding vertex i+1. As a result there are **no duplicate junction points**.
  - total = Σnkpts(1:npaths-1). The last vertex's nk is ignored.
  - The last generated point is overwritten by the final vertex, so the final step is non-uniform (2·dk).
  - `xpoints` accumulates `norm2(dk)` in **fractional** coordinates. There is no reciprocal metric, so distances are not Cartesian.
- Labels: **NOT PRESENT**. Trailing `!X` comments are ignored by list-directed READ.
- `ReadPCPath` uses `NextDataLine`, so it allows comments. `ReadBandsConfig` reads stdin directly.

## 11. OpenMP
- `calc_band_structure`, `calc_band_structure_non_orth` (mod_core.f90 l.116, l.201): `!$omp parallel do default(shared) private(ik,thread_id)` with per-thread pools (`Hk_pool(:,:,tid)`, Sk_pool, w/work/rwork/info pools) sized by `omp_get_max_threads()`.
- `RunUnfolding` (l.999-1119): `!$omp parallel` / `!$omp do` over PC k-points.
  - Per-thread pools: Hk, Sk, Cvec, Cprime, Ssqrt, Usq, Eval, wS, phase.
  - Explicit private scalars. A named critical section `unfold_progress_report`.
  - Solver workspaces are local allocatables inside the called routines, so they are thread-safe.
- Thread id fallback via `!$` sentinels.

## 12. Units / Fermi summary
- Energies are in eV, relative to E_F by default (shift baked into H). S is dimensionless.
- k is fractional everywhere.
- Structure files are in Å and Cartesian (unfold only).
- No unit or E_F metadata is stored. No Bohr appears anywhere in the output path.

## 13. Tests, CI, repository hygiene, examples
- fpm tests: **NOT PRESENT** (auto-tests=false, no `[[test]]`, no test/ directory).
- CI (`.github`): **NOT PRESENT**.
- `.gitignore`: `./build` (a leading `./` pattern is not a valid gitignore path; build/ is in fact ignored by fpm's own `build/.gitignore`), `*.o`, `*.mod`, `*.bin`, `*.zip`.
- Docs: README.md, LICENSE (GPL-3.0-or-later), examples/UnfoldGaNPaperH/VALIDATION.md, `.fortls`.
- Examples (all .hsio files are tiny synthetic ones, 0.6 to 3.3 kB):
  - UnfoldChain1D: Chain_H/S.hsio (760/616 B) and unfold.nml with a manual &ATOMMAP.
  - UnfoldChain1D_Auto: Chain_H/S.hsio (1.1 kB each), pc.nml, sc.nml.
  - UnfoldChain1DShearedM: unfold.nml only (reuses the UnfoldChain1D .hsio files).
  - UnfoldDiatomic: Diatomic_H/S.hsio (1.5 kB/840 B), pc.nml, sc.nml.
  - UnfoldSpinorNonCollinear: Spinor_H/S.hsio (3.3/1.6 kB, atom-major spinor), pc.nml, sc.nml.
  - UnfoldGaNPaperH: OpenMX inputs pc_H.dat/sc_H.dat, UnityCell.nml/SuperCell.nml. The .hsio files referenced by unfold.nml (`openmx/sc_H_H.hsio`, `sc_H_S.hsio`) are NOT PRESENT.
  - Tellurium: Telurio.nml (a "bands" input). Its .hsio files are NOT PRESENT.
  - UnfoldSquare2AtomNonDiag: added in commit 618dc66 but NOT PRESENT in the tree.
- No real (large) DFT .hsio is tracked.

## 14. app/main.f90 / band structure code
- `program LowFoldIn` reads the `&SYSTEM` namelist (ProjectName, Calculation, nonOrtho, HDataFile, SDataFile, OutputDir) plus `&BANDS` from stdin via `config::ReadConfigFromStdin`.
- It calls `core::calc_band_structure` (orthogonal zheev) or `calc_band_structure_non_orth` (zhegv, eigenvalues only).
- Output goes through `bands_io::BandsDataFormat`: `SaveToText` (.bands.dat), `SaveToBinary` (.bands.bin), `SaveToGnuplotBlocks`.
- **Not built as an executable** (it has no `[[executable]]` entry and auto-executables=false). Its modules (core, config, bands_io) are compiled into the library.
- `mod_terminal_ui`: print_box, print_row_split, print_box_close, print_line, UI_WIDTH (UI only).
