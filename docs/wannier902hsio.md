# `wannier902hsio` — Wannier90 Hamiltonian to `.hsio`

`wannier902hsio` is the Wannier90 counterpart of `openmx2hsio`. It keeps
format-specific parsing at the repository boundary and lets the numerical
tools consume the same `.hsio` representation regardless of the originating
electronic-structure code.

```bash
fpm run wannier902hsio -- seedname_hr.dat \
  --out-prefix model \
  --orbs-per-atom "10,10,6,6" \
  --fermi 0.0
```

The command writes `model_H.hsio` and `model_S.hsio`. The latter contains
`S(R) = I delta(R,0)`, because Wannier functions form an orthonormal basis.

Wannier90 degeneracy weights are applied while reading `*_hr.dat`. The stored
blocks therefore follow the LowFoldIn convention directly:

```text
H(k) = sum_R H(R) exp(+i 2 pi k.R)
```

`--orbs-per-atom` is important when atom-resolved tools will be used. Its
positive integers must sum to `num_wann`. Without it, every Wannier function
is represented as one logical atom, which preserves the Hamiltonian but not
the physical atom grouping. For spinor models, include both spin components
in the count for each atom.

`--fermi E_EV` subtracts `E_EV I` from the `R=(0,0,0)` block. It defaults to
zero because the correct reference depends on how the Wannier model was
generated. The input is never overwritten.

The converter does **not** reconstruct real-space Wannier-function amplitudes
from centres, spreads, `*.amn`, `*.mmn`, or `*.chk` files.
