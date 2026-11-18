# `surface_dos` — semi-infinite surface Green function

`surface_dos` computes the retarded Green function of an orthonormal bulk
Hamiltonian represented in `.hsio` format. It is intended as the first
validation stage before real-space LDOS or STM maps are attempted.

The calculation is periodic in the two surface directions and semi-infinite
along `normal_axis`. Unit cells are grouped automatically into a principal
layer thick enough that only neighbouring principal layers couple. The
surface Green function is obtained by iterative Lopez-Sancho decimation.

```fortran
&SURFACE_INPUT
  h_file = 'model_H.hsio'
  output_file = 'model.surface_dos.dat'
  side = 'left'
  normal_axis = 3
  nk1 = 121
  nk2 = 1
  k_origin = 0.0, 0.0, 0.0
  k_span1 = 0.5, 0.0, 0.0
  k_span2 = 0.0, 0.0, 0.0
  energy_min = -0.05
  energy_max = 0.05
  energy_num = 201
  eta = 0.001
  tolerance = 1.0e-10
  max_iter = 100
/
```

Run it with `fpm run surface_dos -- surface_dos.nml`. The columns are the two
grid fractions, energy in eV, total surface DOS in states/eV, and the number
of decimation iterations. Blank lines separate k-points. `side='left'` and
`side='right'` exchange the inter-layer coupling and represent the two
possible terminations of the same bulk ordering.

Current scope and assumptions:

- the semi-infinite model is orthonormal (`S=I`), as produced from Wannier90;
- `eta` is a numerical Lorentzian broadening, not temperature;
- the trace covers one complete principal layer;
- no thermal convolution is applied;
- the output is a surface spectral DOS, not a simulated STM image.

The existing finite-slab `spectral` path retains support for non-trivial
`S(R)` from LCAO calculations. Extending the semi-infinite solver to the
generalized resolvent `z S - H` is a separate step and does not require a
change to `.hsio`.

For a spatial LDOS, the full orbital spectral matrix is available through
`spectral_matrix`. A quantitative Tersoff–Hamann map additionally requires
complex Wannier amplitudes on a real-space grid in the vacuum:

```text
rho(r,E) = sum_(k,alpha,beta)
           phi_(alpha,k)(r) A_(alpha,beta)(k,E) phi*_(beta,k)(r)
```

Wannier centres or orbital labels alone are insufficient for this step.
