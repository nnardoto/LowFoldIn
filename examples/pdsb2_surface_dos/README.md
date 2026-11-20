# PdSb2 surface-DOS validation

This example records the small, reproducible inputs needed to compare
LowFoldIn with an existing WannierTools calculation. It intentionally does
not commit the 282 MB Hamiltonian or any other calculation output.

The source data inspected on `cbagre` are:

```text
/home/cbagre/works/Nathanael/PdSb2/DFT_OpenMX/CWF_Pd-d_Sb-p/PdSb2_hr.dat
/home/cbagre/works/Nathanael/PdSb2/DFT_OpenMX/CWF_Pd-d_Sb-p/
  WT_surface_001/convergence_NP/NP4/wt.in
/home/cbagre/works/Nathanael/PdSb2/DFT_OpenMX/CWF_Pd-d_Sb-p/
  WT_surface_001/convergence_NP/NP4/dos.dat_{l,r,bulk}
```

The model has 88 spinor Wannier functions on four Pd atoms (`5d x 2 spin`)
and eight Sb atoms (`3p x 2 spin`). Convert a local or otherwise accessible
copy without adding it to Git:

```bash
fpm run wannier902hsio -- /path/to/PdSb2_hr.dat \
  --out-prefix PdSb2 \
  --orbs-per-atom "10,10,10,10,6,6,6,6,6,6,6,6"
```

Then run the two terminations:

```bash
fpm run surface_dos -- examples/pdsb2_surface_dos/left_GX.nml
fpm run surface_dos -- examples/pdsb2_surface_dos/right_GX.nml
```

The checked WannierTools reference uses the (001) surface with in-plane
vectors `[100]` and `[010]`, `Eta_Arc=0.001 eV`, 201 energies from -0.05 to
+0.05 eV, 121 points per path segment, and `NP=4`. These namelists reproduce
the first path segment Gamma–X. Compare `PdSb2.left.GX.dat` to the Gamma–X
part of `dos.dat_l`, and analogously for the right side. Absolute intensity
normalisation may differ; energy positions and dispersive features are the
primary validation targets.

No spatial Wannier-function files (`UNK*` or per-WF XSF/CUBE grids) were
found in the inspected project tree. `PdSb2.wout` contains centres and
spreads, while the CUBE files are OpenMX densities/potentials. Consequently,
this example stops at the surface spectral DOS. It must not be described as
a quantitative STM simulation.
