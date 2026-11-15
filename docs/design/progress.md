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
| F1 contracts | done | see git log | contracts.md, test harness |

## Pending / limitations
