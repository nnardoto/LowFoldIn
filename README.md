# LowFoldIn

**LowFoldIn** is a modern Fortran toolkit for post-processing electronic-structure calculations in localized atomic-orbital (LCAO) bases. Its current focus is the calculation of conventional band structures and the unfolding of supercell bands into a primitive-cell Brillouin zone.

The project includes a native OpenMX converter, so an OpenMX `.scfout` file can be transformed directly into the internal Hamiltonian and overlap format used by the band and unfolding solvers.

> **Project status:** research software under active development. The unfolding implementation currently supports diagonal primitive-cell-to-supercell transformation matrices.

## Main features

- Conventional band-structure calculation from real-space Hamiltonian matrices `H(R)`.
- Generalized eigenvalue problem for non-orthogonal bases using `H(k)c = E S(k)c`.
- Band unfolding for LCAO bases with non-orthogonal orbitals.
- Automatic atom mapping between primitive and supercells from atomic positions.
- Optional manual atom mapping for testing and debugging.
- Lorentzian spectral-function output `A(k,E)`.
- OpenMP parallelization over k-points.
- Per-k-point progress and elapsed-time reporting.
- Native conversion of OpenMX 3.9/4.0 `.scfout` files.
- Support for non-spin-polarized, collinear-spin and non-collinear/spinor OpenMX calculations.

## Scientific basis

The unfolding implementation follows the formalism described in:

> J. Quan, N. Rybin, M. Scheffler and C. Carbogno,  
> *Efficient band structure unfolding with atom-centered orbitals: General theory and application*,  
> **Physical Review B 113**, 085112 (2026).  
> DOI: https://doi.org/10.1103/7xym-7388

## Requirements

- A Fortran compiler with Fortran 2023 support, such as recent `gfortran`.
- [Fortran Package Manager (fpm)](https://fpm.fortran-lang.org/), version 0.13 or newer recommended.
- BLAS.
- LAPACK.
- OpenMP support for parallel execution.

On Arch Linux, the basic dependencies can be installed with:

```bash
sudo pacman -S gcc-fortran blas lapack
```

Install `fpm` using your preferred package source or the instructions from the official fpm documentation.

## Building

Clone the repository and enter its directory:

```bash
git clone https://github.com/SEU-USUARIO/LowFoldIn.git
cd LowFoldIn
```

Build the project:

```bash
fpm build --profile release
```

The `fpm.toml` file defines optimization, diagnostics and OpenMP flags through the `meus_flags` feature. Depending on the installed fpm version, it may be necessary to activate the feature explicitly or pass the compiler flags directly:

```bash
fpm build --profile release --flag "-std=f2023 -O3 -Wall -Wextra -fopenmp"
```

For runtime debugging, a useful alternative is:

```bash
fpm build --flag "-std=f2023 -O0 -g -Wall -Wextra -fcheck=all -fbacktrace -fopenmp"
```

## Executables

The package builds three command-line programs.

| Executable | Purpose |
|---|---|
| `LowFoldIn` | Conventional band-structure calculation. Reads its input from standard input. |
| `openmx2hsio` | Converts an OpenMX `.scfout` file into the internal `.hsio` format. |
| `unfold` | Unfolds a supercell band structure into the primitive-cell Brillouin zone. |

They can be launched through fpm:

```bash
fpm run LowFoldIn < bands.nml
fpm run openmx2hsio -- calculation.scfout
fpm run unfold -- unfold.nml
```

After building, the corresponding binaries are also available inside the fpm build directory.

---

## Workflow with OpenMX

A typical workflow is:

1. Run an OpenMX calculation and retain the `.scfout` file.
2. Convert the Hamiltonian and overlap matrices with `openmx2hsio`.
3. Prepare primitive- and supercell structure files.
4. Define the primitive-cell k-path in `unfold.nml`.
5. Run `unfold`.
6. Plot either the discrete spectral weights or the broadened spectral function.

### 1. Convert an OpenMX `.scfout` file

Basic usage:

```bash
fpm run openmx2hsio -- system.scfout
```

By default, the output prefix is derived from the `.scfout` filename.

Available options:

```text
openmx2hsio <file.scfout> [--out-prefix PREFIX] \
             [--fermi E_F_HARTREE] [--no-fermi-shift]
```

Example:

```bash
fpm run openmx2hsio -- supercell.scfout --out-prefix data/supercell
```

For a non-spin-polarized calculation, the converter generates:

```text
data/supercell_H.hsio
data/supercell_S.hsio
```

For collinear spin, it generates:

```text
data/supercell_H_up.hsio
data/supercell_H_dn.hsio
data/supercell_S.hsio
```

Run the unfolding separately for the up and down Hamiltonians, using the same overlap matrix and atom mapping.

For a non-collinear/spinor calculation, the converter writes atom-major spinor matrices compatible with the unfolding code.

By default, the Hamiltonian is shifted so that the Fermi energy is zero. The Fermi level is read from `ChemP` in the `.scfout`, unless overridden with `--fermi`. Use `--no-fermi-shift` to preserve the original energy reference.

---

## Band-structure calculation

The main `LowFoldIn` executable reads a `&SYSTEM` namelist followed by a `&BANDS` block from standard input.

Example `bands.nml`:

```fortran
&SYSTEM
    ProjectName = "example"
    Calculation = "bands"
    nonOrtho    = .true.
    HDataFile   = "data/example_H.hsio"
    SDataFile   = "data/example_S.hsio"
    OutputDir   = "results"
/

&BANDS
4
0.0  0.0  0.0   50
0.5  0.0  0.0   50
0.5  0.5  0.0   50
0.0  0.0  0.0    1
/
```

Run it with:

```bash
mkdir -p results
fpm run LowFoldIn < bands.nml
```

The program produces:

```text
results/example.bands.dat
results/example.bands.bin
results/example.bands.gnuplot.dat
```

- `.bands.dat`: text output with path coordinate and band energies.
- `.bands.bin`: binary output intended for subsequent processing.
- `.bands.gnuplot.dat`: bands separated into blocks for direct use with Gnuplot.

Set `nonOrtho = .true.` whenever an overlap matrix `S(R)` is required. In this mode, the code solves the generalized Hermitian eigenvalue problem.

---

## Band unfolding

The dedicated `unfold` executable reads an input file containing:

1. `&UNFOLD` configuration;
2. an optional manual `&ATOMMAP` block;
3. the primitive-cell path in `&PCPATH`.

Automatic atom mapping from structure files is recommended.

### Minimal unfolding input

```fortran
&UNFOLD
    ProjectName = "example"
    HDataFile   = "data/supercell_H.hsio"
    SDataFile   = "data/supercell_S.hsio"
    OutputDir   = "results"

    M(1,:) = 4, 0, 0
    M(2,:) = 0, 1, 0
    M(3,:) = 0, 0, 1

    PCStructureFile = "structures/primitive.nml"
    SCStructureFile = "structures/supercell.nml"
    MatchTolerance  = 1.0e-3

    WriteSpectralFunction = .true.
    Broadening = 0.02
    Emin = -5.0
    Emax =  5.0
    nE   = 1000
/

&PCPATH
4
0.0  0.0  0.0   80
0.5  0.0  0.0   80
0.5  0.5  0.0   80
0.0  0.0  0.0    1
/
```

Run:

```bash
mkdir -p results
export OMP_NUM_THREADS=8
fpm run unfold -- unfold.nml
```

During diagonalization, a message is printed immediately after each k-point is completed:

```text
k-point 7 completed | progress 5 / 120 (4.2%) | time: 12.384 s
```

Because the calculation is parallel, completed k-points may be reported out of numerical order. The progress counter always represents the total number of finished points.

### Transformation matrix `M`

The matrix relates primitive-cell and supercell lattice vectors. It must be written explicitly by rows:

```fortran
M(1,:) = 4, 0, 0
M(2,:) = 0, 2, 0
M(3,:) = 0, 0, 1
```

The current implementation only accepts diagonal matrices:

```text
M = diag(m1,m2,m3)
```

Non-diagonal matrices are rejected with an explicit error instead of producing a potentially incorrect result.

### Structure-file format

The primitive- and supercell structure files use the same simple text format.

```fortran
&STRUCTURE
    LatticeVectors(1,:) = 4.4566,  0.0000, 0.0000
    LatticeVectors(2,:) = -2.2283, 3.8596, 0.0000
    LatticeVectors(3,:) = 0.0000,  0.0000, 5.9280
/

&SPECIES
1
Te  13
/

&ATOMS
3
Te  1.1732  0.0000  1.9760
Te  3.2834  1.9270  3.9520
Te  2.2283  3.8596  0.0000
/
```

Important details:

- Lattice vectors are given in ångström.
- Each **row** of `LatticeVectors` is one lattice vector, following the convention commonly used by POSCAR, ASE and CIF tools.
- Atomic coordinates are Cartesian and given in ångström.
- Atoms in the supercell structure file must appear in the same order as in the `.hsio` file.
- The `&SPECIES` block is optional, but strongly recommended because it allows the program to verify the number of orbitals per atom against the actual `.hsio` data.
- Fortran namelists use parentheses, such as `M(1,:)`, not square brackets.
- Do not flatten a non-orthogonal lattice matrix into a single assignment. Fortran fills arrays in column-major order, which may silently transpose the intended lattice.

### Automatic atom mapping

When both `PCStructureFile` and `SCStructureFile` are present, the program constructs the map

```text
supercell atom -> primitive-cell atom + primitive-cell translation
```

using:

- atomic species;
- atomic positions;
- periodic minimum-image distances;
- the transformation matrix `M`;
- optionally, the number of orbitals declared in `&SPECIES`.

`MatchTolerance` is expressed in primitive-cell fractional coordinates.

### Manual atom mapping

For small systems or debugging, omit the two structure filenames and provide `&ATOMMAP` directly:

```fortran
&ATOMMAP
1  0 0 0
2  0 0 0
1  1 0 0
2  1 0 0
/
```

Each line contains:

```text
primitive_atom_index  n1  n2  n3
```

The lines must follow exactly the atom order stored in the `.hsio` file. Translation indices are zero-based and must satisfy `0 <= ni < M(i,i)`.

### Unfolding outputs

The program writes:

```text
results/example.unfold.dat
results/example.unfold_spectral.dat
```

`example.unfold.dat` contains triples:

```text
x   energy   weight
```

It is suitable for a scatter plot in which marker intensity, size or color represents the unfolded spectral weight.

`example.unfold_spectral.dat` contains the Lorentzian-broadened spectral function `A(k,E)` on the energy grid defined by `Emin`, `Emax`, `nE` and `Broadening`.

The code also reports the unfolding sum rule. For `p` primitive-cell orbitals, the expected value is:

```text
sum_N W = p
```

The largest observed deviation is printed after the calculation. A large deviation generally indicates an inconsistent transformation matrix, atom mapping, structure ordering or orbital count.

---

## Build profiles

LowFoldIn requires `fpm` 0.13 or newer and is configured for GNU Fortran 14 or newer. The compiler version is not automatically enforced by the manifest, so select it explicitly when needed.

```bash
# Full runtime checks and debugging information
fpm build --profile debug --compiler gfortran-14

# Balanced optimization
fpm build --profile balanced --compiler gfortran-14

# Portable production build
fpm build --profile release --compiler gfortran-14

# Optimized for the machine performing the compilation
fpm build --profile native --compiler gfortran-14

# AMD Ryzen 9 7950X / Zen 4
fpm build --profile ryzen7950x --compiler gfortran-14
```

The `ryzen7950x-fast` profile additionally enables `-ffast-math` and loop unrolling. It may change floating-point behavior and should only be used after validating eigenvalues and spectral weights against the conservative `ryzen7950x` profile.

## Parallel execution

The expensive k-point loops are parallelized with OpenMP. Set the thread count before running:

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
fpm run unfold -- unfold.nml
```

For reproducibility or debugging, run serially:

```bash
OMP_NUM_THREADS=1 fpm run unfold -- unfold.nml
```

Each thread owns its working matrices for `H(k)`, `S(k)`, eigenvectors and auxiliary data. Results from different k-points are independent.

## Repository structure

```text
LowFoldIn/
├── app/
│   └── main.f90                 # Conventional band-structure executable
├── src/
│   ├── mod_bands_io.f90         # Band-result writers
│   ├── mod_config.f90           # Main-program input parser
│   ├── mod_core.f90             # Fourier transforms and eigensolvers
│   ├── mod_hs_io.f90            # Internal H(R)/S(R) binary format
│   ├── mod_openmx_convert.f90    # OpenMX matrix assembly and conversion
│   ├── mod_scfout_io.f90        # OpenMX .scfout reader
│   ├── mod_terminal_ui.f90       # Terminal formatting
│   ├── mod_types.f90             # Numerical kinds and shared types
│   └── mod_unfold.f90            # Unfolding implementation
├── tools/
│   ├── openmx2hsio.f90          # OpenMX conversion CLI
│   └── unfold.f90               # Unfolding CLI
├── fpm.toml
├── LICENSE
└── README.md
```

## Current limitations

- Only OpenMX `.scfout` format version 3 is supported, corresponding to OpenMX 3.9/4.0.
- The unfolding transformation matrix `M` must currently be diagonal.
- The structure-file parser expects Cartesian coordinates in ångström.
- The supercell atom order must match the order represented in the converted `.hsio` file.
- There is not yet a stable public library API; the command-line programs are the primary interface.
- Automated regression tests and packaged examples should be added before treating the project as production software.

## Recommended additions before a public release

The following repository additions would make the project easier to reproduce and cite:

- `CITATION.cff`;
- small test datasets under `examples/`;
- plotting scripts for `.bands.gnuplot.dat`, `.unfold.dat` and `.unfold_spectral.dat`;
- continuous integration for serial and OpenMP builds;
- numerical regression tests for the unfolding sum rule and spin channels.

## Citation

When using the unfolding functionality, cite the theoretical reference above. A dedicated citation for LowFoldIn can be added when the repository receives an archived release and DOI.

## License

LowFoldIn is free software distributed under the **GNU General Public License, version 3 or any later version** (`GPL-3.0-or-later`). See [`LICENSE`](LICENSE) for the complete license text.

## Author

Nathanael Nardoto Batista

## Contributing

Issues and pull requests are welcome. For numerical changes, please include:

- the compiler and compiler version;
- BLAS/LAPACK implementation;
- OpenMP thread count;
- a minimal input reproducing the result;
- the expected and obtained sum-rule deviations.
