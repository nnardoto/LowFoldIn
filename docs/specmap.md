# `specmap` — broadened spectral maps from `.spectral.bin`

`specmap` is the **Represent** step of the surface spectral pipeline
(`atomselect` → `spectral` → `specmap` → gnuplot). It reads the discrete
energies and projection weights written by [`spectral`](spectral.md) and
broadens them into continuous maps A_P(k,E) on a (k, E) grid, for single
regions or ± combinations of regions, plus a gnuplot include with the k-path
ticks.

`specmap` reads **only** the `.spectral.bin` file: no HSIO, no structure, no
LAPACK, no diagonalisation. Changing the kernel, width, energy window, maps
or occupation only requires re-running `specmap`.

During execution it reports the input-validation, map-generation and metadata
stages. Each requested map shows its position in the queue and elapsed time;
the final panel reports the number of maps and total wall time.

## Usage

```bash
specmap <specmap.nml>
fpm run specmap -- specmap.nml
```

Exactly one argument. On error the message is printed prefixed
`[specmap error]` and the exit code is 1. All map expressions are validated
before any file is written.

## Input: `&SPECMAP`

| field | type | default | meaning |
|---|---|---|---|
| `spectral_file` | string | `''` | the `.spectral.bin` file (**required**) |
| `spin_channel` | integer | `1` | 1..nspin (collinear: 1 = up, 2 = down) |
| `kernel` | string | `'lorentzian'` | `'lorentzian'` or `'gaussian'` (case-insensitive) |
| `width_ev` | real (eV) | `0.05` | Lorentzian: HWHM η; Gaussian: standard deviation σ; must be > 0 |
| `emin`, `emax` | real (eV) | `-3.0`, `3.0` | energy window, `emax > emin`, in the file's energy reference |
| `ne` | integer | `601` | number of energy points, ≥ 2 |
| `maps` | strings | — | up to 32 map expressions (at least one) |
| `occupation` | logical | `.false.` | multiply by the Fermi–Dirac occupation |
| `temperature_k` | real (K) | `0.0` | temperature for the occupation, ≥ 0 |
| `mu_ev` | real (eV) | `0.0` | chemical potential, used only if `mu_set` |
| `mu_set` | logical | `.false.` | `mu_ev` is given explicitly |
| `output_prefix` | string | `'specmap'` | output file prefix (must not be empty) |

Energy grid: `E_i = emin + (i−1)(emax−emin)/(ne−1)`, `i = 1..ne`.

### Map expressions

A map expression is a list of region names (as stored in the spectral file,
including `total`) joined by `+` or `-`; a leading `-` is allowed and spaces
are ignored:

```fortran
maps = 'total', 'top', 'top-bottom', 'top+bottom', 'a+b-c'
```

- The map is the linear combination with coefficients ±1 of the region
  weights (combined before broadening, which is identical by linearity).
  There is **no normalisation** and the sign is preserved, so difference
  maps such as `top-bottom` are positive where `top` dominates and negative
  where `bottom` dominates.
- A repeated region adds its coefficients (`top+top` = 2·top).
- Errors: empty expression, empty term (`top+`, `++`, `+top`), unknown region
  (the message lists the available names).
- Summing **overlapping** regions double-counts the shared atoms.

### Occupation and chemical potential

With `occupation = .true.` the intensity is `I(k,E) = f(E)·A(k,E)` with

```text
f(E) = 1 / (1 + exp((E − μ)/(k_B T))),   k_B = 8.617333262e-5 eV/K
```

evaluated in an overflow-safe way; `T = 0` gives a step with `f(μ) = 1/2`.
The value of μ is chosen as follows (E_F is **never** used implicitly):

| file energy reference | `mu_set = .true.` | `mu_set = .false.` |
|---|---|---|
| `fermi` | μ = `mu_ev` | μ = 0 (= E_F by definition) |
| `absolute` | μ = `mu_ev` | error |
| `unknown` | μ = `mu_ev` | error |

## Equations

```text
A_P(k,E) = Σ_n w_n^P(k) K(E − ε_n(k))
```

| kernel | `width_ev` means | K(x) |
|---|---|---|
| Lorentzian | half width at half maximum η | (1/π) η / (x² + η²) |
| Gaussian | standard deviation σ (**not** HWHM; FWHM = 2√(2 ln 2) σ ≈ 2.355 σ) | exp(−x²/2σ²) / (σ√(2π)) |

Both kernels integrate to 1, so a map integrated over all energies gives the
summed weights. Weights are used with their sign (Mulliken weights may be
slightly negative). Lorentzian tails are always computed in full; Gaussian
contributions beyond 12σ are skipped (negligible in double precision).
Because the grid is finite, intensity outside `[emin, emax]` is lost.

Units: k-distance in 1/Å, energy in eV, intensity in states/eV (weighted by
Mulliken or Löwdin populations, whichever `spectral` used).

## Output

### `<prefix>.<tag>.map.dat`

One file per map expression. `tag` is the expression with spaces removed, a
leading `-` replaced by `minus_`, then `+` → `_plus_` and `-` → `_minus_`
(e.g. `top-bottom` → `model.top_minus_bottom.map.dat`).

Header (lines starting with `#`):

```text
# specmap spectral map (broadened projection weights; NOT a photoemission simulation)
# source file: model.spectral.bin
# source (from file): H=model_H.hsio S=model_S.hsio
# map expression: top   spin channel: 1
# method: mulliken
# kernel: lorentzian, width = 0.05 eV (HWHM eta)
# energy window: [-3.5, 3.5] eV, ne = 701
# occupation: none
# units: k_distance 1/Angstrom, energy eV, intensity states/eV (Mulliken/Lowdin-weighted)
# energy reference: fermi (E - E_F), E_F = 0.0 eV
# regions combined with coefficients +/-1, no normalisation; overlapping regions in a sum may double-count atoms
# columns: k_distance energy_eV intensity
```

Then, for every k-point, `ne` rows `k_distance energy_eV intensity`
followed by one blank line (gnuplot `pm3d` block format; `nk·ne` data lines
and `nk` blank lines).

### `<prefix>.kticks.gp`

A gnuplot include for the k axis and energy reference, e.g.

```gnuplot
set xrange [0:2.52815853E+000]
set xtics ("X" 0.00000000E+000, "{/Symbol G}" 1.04719755E+000, "M" 2.52815853E+000)
if (!exists("VLINE_COLOR")) VLINE_COLOR = "white"
set arrow from 1.04719755E+000, graph 0 to 1.04719755E+000, graph 1 nohead dt 2 lw 1 lc rgb VLINE_COLOR front
EF = 0.0
ELABEL = "E - E_F (eV)"
```

- Vertex labels at their `kdist`; `G`/`Gamma` become `{/Symbol G}`.
- Dashed vertical lines at interior vertices (colour `VLINE_COLOR`, default
  white; define it before loading to change it).
- `EF` = 0 for the `fermi` reference, E_F for `absolute` with a known E_F,
  otherwise `NaN`. `ELABEL` = `"E - E_F (eV)"`, `"E (eV)"` or
  `"E (eV, reference unknown)"`.

## Plotting with gnuplot

Generic scripts in `tools/gnuplot/` (each writes `<out>.png` and
`<out>.pdf`):

```bash
gnuplot -e "datafile='slab.top.map.dat'; ticks='slab.kticks.gp'; out='slab_top'" tools/gnuplot/map_positive.gp
gnuplot -e "datafile='slab.top_minus_bottom.map.dat'; ticks='slab.kticks.gp'; out='slab_diff'" tools/gnuplot/map_difference.gp
gnuplot -e "datafile1='slab.top.map.dat'; datafile2='slab.bottom.map.dat'; ticks='slab.kticks.gp'; out='slab_cmp'" tools/gnuplot/map_compare.gp
```

- `map_positive.gp`: sequential palette (non-negative maps); optional
  `title`, `cbmax`.
- `map_difference.gp`: divergent blue–white–red palette centred at 0,
  symmetric range `[-m:m]`; optional `title`, `cbmax`.
- `map_compare.gp`: two maps side by side with a common colour range;
  optional `title1`, `title2`, `cbmax`.

## Approximations and limitations

- The maps are **broadened projection weights**, not photoemission
  intensities: no matrix elements, photon energy/polarisation, final
  states, k_z broadening, experimental resolution or many-body
  self-energies. The broadening is an artificial, user-chosen width.
- Only k-path (line) maps; no constant-energy kx × ky maps.
- No normalisation of combined maps; overlapping regions double-count.
- Region names containing `+` or `-` could not be addressed (such names are
  rejected by `atomselect`/`regions.nml` anyway).
- `spin_channel` selects one channel; up+down sums of a collinear file are
  not produced.

## Band structure, weights, spectral function, maps and ARPES

See the corresponding section of [spectral.md](spectral.md#band-structure-weights-spectral-function-maps-and-arpes):
bands are ε_n(k); weights w_n^P(k) are populations of region P; A_P(k,E) is
those weights broadened by an artificial kernel; an intensity map is A_P (or a
± combination, optionally × f(E)) on a grid; experimental ARPES additionally
involves matrix elements, final states, resolution and many-body effects,
none of which are modelled.

## Tests

`fpm test test_specmap` builds a synthetic spectral data set in code and
checks: kernel normalisation, intensity conservation, re-broadening of the
same file with a different width, map combinations (`+`, `-`, repeated and
unknown regions), negative weights kept with their sign, Fermi–Dirac
(including T = 0 and overflow safety), the map file layout, file-name tags,
the kticks file and rejection of a truncated spectral file.

## Worked example

[`examples/surface_spectral`](../examples/surface_spectral/README.md) runs two
`specmap` inputs on the same `model.spectral.bin`:

```fortran
&SPECMAP                                   ! specmap.nml
    spectral_file = 'model.spectral.bin'
    spin_channel  = 1
    kernel        = 'lorentzian'
    width_ev      = 0.05
    emin = -3.5, emax = 3.5, ne = 701
    maps          = 'total', 'top', 'bottom', 'interior', 'top-bottom', 'top+bottom'
    occupation    = .false.
    output_prefix = 'model'
/
```

and `specmap_gauss.nml` (Gaussian, σ = 0.10 eV, `output_prefix = 'model_gauss'`),
then plots them with `plot_spectral.gp`.
