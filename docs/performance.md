# Performance notes — `spectral`

## Building with OpenMP

```sh
fpm build --profile release        # -O3 -fopenmp (see fpm.toml [profiles])
fpm run spectral --profile release -- spectral.nml
```

A plain `fpm build` has no `-fopenmp`: the `!$` lines are ignored and the
k-loop runs serially with identical results.

## Environment (examples, not requirements)

```sh
export OMP_NUM_THREADS=4
# Avoid oversubscription: one OpenMP thread per k-point already, so keep the
# BLAS/LAPACK library single-threaded inside each thread, e.g.
export OPENBLAS_NUM_THREADS=1      # OpenBLAS
export VECLIB_MAXIMUM_THREADS=1    # Apple Accelerate / vecLib
```

The program never sets environment variables; it prints the number of
threads available at start-up and repeats this hint when it is > 1.

## What is parallelized

- `run_spectral`: the loop over the flat index `(ispin, ik)` is an
  `!$omp do schedule(dynamic,1)`. Each iteration builds H(k), S(k), solves
  the generalized eigenproblem and computes the projection weights
  (`process_kpoint`), then writes only `sd%energy(:,ik,ispin)` and
  `sd%weight(:,ik,:,ispin)`. Results are therefore independent of the
  thread count (checked bitwise in `test_spectral_engine`, 1 vs 4 threads).
- Errors inside the loop are stored per iteration; after the loop the one
  with the lowest index is reported (deterministic). Normalisation /
  negative-weight warnings are printed after the loop (`check_weights`).
- Reading inputs, the k-path and writing the output stay serial.

## Memory

Per thread: Hk, Sk, C, SC (4 n² complex) + E (n reals) + W (n·ngroup reals),
allocated once per thread. The LAPACK wrappers in `spectral_projection`
additionally allocate their own temporaries per call (copies of H and S for
`zhegv`, 2–5 n² complex for Löwdin). Eigenvectors are never stored for all
k-points: only energies (n·nk·nspin) and weights (n·nk·ngroup·nspin) are kept.

## Timing report

All values come from `system_clock` (wall clock). `read`, `k-loop`, `write`
and `total` are wall times. `build`, `diagonalization` and `projection` are
**accumulated** over threads (sum of each thread's wall time in that phase,
thread-seconds); with several threads they can exceed the k-loop wall time.

## Measured timings

Machine: Apple M2 (8 cores), macOS 27.0.1, GNU Fortran 16.1.0, fpm 0.13.0,
`--profile release`, LAPACK/BLAS = Apple Accelerate,
`VECLIB_MAXIMUM_THREADS=1`. Case: `examples/UnfoldSpinorNonCollinear`
(n = 18 orbitals, non-collinear, Mulliken, 1 region), straight path with
nk = 2000. Three runs each; values in seconds.

| OMP_NUM_THREADS | k-loop wall | build (acc.) | diag (acc.) | proj (acc.) | total wall |
|---|---|---|---|---|---|
| 1 | 0.1117 / 0.0975 / 0.0978 | 0.0045 / 0.0034 / 0.0035 | 0.0976 / 0.0864 / 0.0866 | 0.0094 / 0.0075 / 0.0075 | 0.1137 / 0.1004 / 0.0990 |
| 4 | 0.0329 / 0.0328 / 0.0334 | 0.0046 / 0.0045 / 0.0046 | 0.1146 / 0.1141 / 0.1162 | 0.0117 / 0.0115 / 0.0120 | 0.0342 / 0.0340 / 0.0347 |

The `.spectral.bin` files from the 1- and 4-thread runs are byte-identical.

**Caveat:** the system is tiny (18×18 matrices, ~50 µs per k-point), so these
numbers are not representative of production slabs, where the O(n³)
diagonalization dominates and memory bandwidth / BLAS threading matter.
