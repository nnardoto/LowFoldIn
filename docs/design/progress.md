# Surface spectral analysis — progress log

Branch: `feat/surface-spectral-analysis` (base: `main` @ 692012e)

## Decisions
- C0.1: HSIO has no geometry → reuse the existing structure `.nml` file
  (&STRUCTURE/&SPECIES/&ATOMS). HSIO binary format untouched.
- Energy: `openmx2hsio` shifts H by −E_F·S by default and does not store E_F.
  `spectral` never shifts; the user declares `energy_reference`.
- Spin: no flag in HSIO; `spin_mode` declared in `spectral.nml`. Noncollinear
  ordering is atom-major [up | down] per atom; regions take all Orbs(a).
- k-path: new `kpath` module (reciprocal metric, labels); `unfold::ReadPCPath`
  untouched to keep `unfold` output unchanged.
- Solver: new non-stopping `generalized_eigensolve` (zhegv on copies);
  `unfold` private solver left unchanged.
- Plain `fpm build` has no OpenMP; profiles add `-fopenmp`.

## Fichas
| Ficha | Status | Commit | Notes |
|---|---|---|---|
| F0 survey | done | (F1 commit) | docs/design/architecture-survey.md |
| F1 contracts | done | dafe60b | contracts.md, test harness |
| F2a structure_io | done | 4660526 | 42 checks |
| F2b regions/selection | done | eee8501 | 70 checks |
| F3 spectral_data | done | b12e32f | 82 checks |
| F5 atomselect CLI | done | 6914a99 | 37 CLI checks (test/test_atomselect_cli.sh) |
| F4 projection + kpath | done | 470061a | 45 + 24 checks |
| F6 spectral exe | done | d844f20 | 131 checks; Diatomic manual run max dev 2e-15 |
| F9 Löwdin | done | 4859668 | |
| F7 specmap (+F10 kernels/ops/FD) | done | 51a0cc4 | 59 checks, gnuplot scripts rendered |
| F8 geometry/layers | done | (see log) | 45 checks + CLI 92 |
| F11 OpenMP | done | (see log) | bitwise identical 1 vs 4 threads; docs/performance.md |

## Pending / limitations
- macOS/Accelerate: `--profile debug` (-ffpe-trap) crashes inside libBLAS
  (SIGILL) in test_projection; not reproduced without the trap flag.
- Species labels > 8 chars truncated in regions.nml (documented in contracts §2).
- Corrupt HSIO → existing ReadFile stops (no ierr path).
- Projection wrappers allocate O(n^2) temporaries per call (negligible vs
  O(n^3) diagonalisation); caller-owned workspaces are a possible refinement.
- Diatomic example has S close to I, so Mulliken vs Löwdin barely differ there.
