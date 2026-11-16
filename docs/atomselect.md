# `atomselect` — named atomic regions

`atomselect` is the **Select** step of the surface spectral pipeline
(`atomselect` → `spectral` → `specmap`). It reads a structure file, reports
atoms and layers, and creates or extends a `regions.nml` file holding named
lists of atom indices (for example `top`, `bottom`, `interior`). The regions
are later used by [`spectral`](spectral.md) to compute per-region projection
weights.

`atomselect` never reads HSIO files and never computes anything electronic:
it only works with atom indices and Cartesian positions.

## Structure file

The structure file is the format already used by `unfold` (rows of
`LatticeVectors` are a1, a2, a3 in Å; atoms in Cartesian Å):

```fortran
&STRUCTURE
    LatticeVectors(1,:) = 3.0, 0.0, 0.0
    LatticeVectors(2,:) = 0.0, 3.0, 0.0
    LatticeVectors(3,:) = 0.0, 0.0, 40.0
/
&SPECIES            ! optional
2
A 2                 ! label, nOrb counted in the HSIO basis (includes spinor doubling)
B 1
/
&ATOMS
8
A 0.0 0.0 13.0      ! label x y z (Cartesian Å), same order as the HSIO
B 0.0 0.0 15.0
...
/
```

- Atom indices are 1-based and follow the order of the file, which must be
  the HSIO atom order. Atoms are never reordered.
- The per-atom label is used as given (no element symbols are invented).
  Labels are stored with at most 8 characters.
- Lines starting with `!` or `#` and blank lines are skipped inside
  `&SPECIES` / `&ATOMS`.
- `&SPECIES` is not used by `atomselect`; `spectral` uses it (when present)
  to check `nOrb(label(a)) == Orbs(a)` against the HSIO.

## Command-line synopsis

```text
atomselect inspect     <structure.nml> [--xyz out.xyz] [--table out.txt]
atomselect fingerprint <structure.nml>
atomselect select      <structure.nml> --regions <regions.nml> --name NAME --atoms "LIST"
atomselect select      <structure.nml> --regions <regions.nml> --name NAME --normal nx ny nz [--smin A] [--smax B]
atomselect select      <structure.nml> --regions <regions.nml> --name NAME --box "xlo xhi ylo yhi zlo zhi"
atomselect layers      <structure.nml> --normal nx ny nz [--tol T] [--refs "c1,c2,..."] [--table out] [--regions file]
atomselect surface     <structure.nml> --regions <file> --normal nx ny nz [--tol T] [--refs "c1,..."]
                       [--top N] [--bottom N] [--top-name top] [--bottom-name bottom]
atomselect validate    <structure.nml> <regions.nml> [--allow-no-fingerprint]
atomselect help        (also -h, --help; running without arguments prints the same text)
```

The structure file (and, for `validate`, the regions file) are positional
arguments right after the subcommand; options may follow in any order.
Exit code is 0 on success and 1 on any error; error messages are printed on
standard error prefixed with `[atomselect error]`. Warnings are prefixed
with `[atomselect warning]` and go to standard output.

### `inspect`

Prints a header (file name, `nat`, the three lattice vectors in Å, the
fingerprint) and the table

```text
 Index  Element           X           Y           Z
```

in the original atom order (Cartesian Å).

| option | meaning |
|---|---|
| `--table out.txt` | also write the header + table to a file |
| `--xyz out.xyz` | also write an XYZ file (labels as given, Cartesian Å) |

### `fingerprint`

Prints the 16-hex-character structure fingerprint (see below) and nothing
else. Useful in scripts.

### `select`

Adds **one** region named `--name` to the file given by `--regions`.
Exactly one selection mode must be used:

| mode | options | selects |
|---|---|---|
| numeric | `--atoms "LIST"` | the listed indices |
| projection range | `--normal nx ny nz` with `--smin A` and/or `--smax B` | atoms with `A ≤ s ≤ B`, where `s = r·n̂` (Å) |
| Cartesian box | `--box "xlo xhi ylo yhi zlo zhi"` | atoms with `lo ≤ x,y,z ≤ hi` (Cartesian Å, inclusive) |

- `LIST` syntax: comma-separated items, each an index or an inclusive range
  `a-b` (`a ≤ b`), spaces optional: `"1-12"`, `"15, 17, 22"`,
  `"1-10, 15-20"`. Errors: empty list, non-integer token, reversed range,
  index < 1 or > nat, duplicates (reported, never silently removed).
- `--normal` takes three Cartesian components (three separate arguments);
  it is normalised internally (zero length is an error). Omitting `--smin`
  or `--smax` leaves that side unbounded, but at least one of them is
  required; `--smin > --smax` is an error. `--smin`/`--smax` without
  `--normal` is an error.
- `--box` takes six limits in one quoted argument; `*` means unbounded on
  that side. A lower limit greater than the upper limit is an error.
- An empty selection is an error; `--atoms`, `--normal` and `--box` are
  mutually exclusive.

The atom list is stored sorted ascending. The program prints, e.g.,
`region 'top': 1 atoms [8] added to regions.nml (created)`.

### `layers`

Projects every atom on the normal, `s_i = r_i · n̂` (Å), and groups atoms
into layers. Prints the header, the normal (normalised), the tolerance, a
table

```text
 Index  Element           X           Y           Z   Projected   Layer
```

and a summary (`# nlayer`, then for each layer its centre = mean `s`, atom
count and atom list).

| option | default | meaning |
|---|---|---|
| `--normal nx ny nz` | required | Cartesian direction (normalised internally) |
| `--tol T` | 0.5 Å | layer tolerance, must be > 0 |
| `--refs "c1,c2,..."` | — | manual layer centres (Å, along n̂): each atom is assigned to the nearest centre; error if an atom is farther than `T` from every centre, or if centres are duplicated |
| `--table out` | — | also write the table to a file |
| `--regions file` | — | also add one region per layer, named `layer_001`, `layer_002`, … |

Layer detection is single-linkage along `s`: atoms are sorted by `s` and a
new layer starts whenever the gap to the previous atom exceeds `T`. Layers
are numbered by **ascending `s`** (layer 1 is the lowest along +n̂). With
`--refs` the layers are numbered by ascending centre, so the numbering is
the same.

### `surface`

Detects layers exactly as `layers` (same `--normal`, `--tol`, `--refs`) and
adds the outermost layers as regions:

| option | default | meaning |
|---|---|---|
| `--regions file` | required | regions file to create/extend |
| `--normal nx ny nz` | required | Cartesian direction |
| `--top N` | — | region from the `N` layers with the highest `s` |
| `--bottom N` | — | region from the `N` layers with the lowest `s` |
| `--top-name NAME` | `top` | name of the top region (requires `--top`) |
| `--bottom-name NAME` | `bottom` | name of the bottom region (requires `--bottom`) |
| `--tol T` | 0.5 Å | as for `layers` |
| `--refs "..."` | — | as for `layers` |

At least one of `--top`/`--bottom` is required; each must be ≥ 1 and
`top + bottom > nlayer` is an error. "Top" always means highest `s` along
**+n̂**: with `--normal 0 0 -1` the "top" region is the one with the lowest
`z`.

### `validate`

Reads the regions file, validates it against the structure and prints:
`nat`, fingerprint status, every region (name, number of atoms, compressed
list), pairwise overlaps (`overlap: 'a' and 'b' share N atoms`), the number
of atoms not covered by any region, whether the regions form a partition
(`partition   : yes|no`, i.e. every atom in exactly one region) and
`valid       : yes`.

`--allow-no-fingerprint` accepts a regions file whose `fingerprint` is empty
(a warning is printed). Without it such a file is rejected.

## Normal direction and non-orthogonal cells

- The surface normal is **never assumed to be z**. `--normal` is a Cartesian
  vector; any direction (e.g. `1 1 0`) can be used, and the tests include
  tilted slabs and `-z` normals.
- Positions are projected in Cartesian space, so non-orthogonal cells need
  no special treatment for the projection itself.
- **Periodic split warning.** Projected coordinates are not wrapped into the
  cell. Before detecting layers `atomselect` checks, for every lattice
  vector `a_i` with `|a_i·n̂| > T` (in-plane vectors are skipped), whether
  the periodic image of the lowest atom would join the top of the slab. If
  so it prints

  ```text
  [atomselect warning] the slab appears to cross the periodic cell boundary along the normal ...
  ```

  and continues **without re-wrapping anything**: the detected layers may
  then be split. Translate the structure so the slab is not cut by the
  boundary, or give explicit centres with `--refs`.

## `regions.nml`

```fortran
&REGION_FILE
    format_version = 1
    natoms         = 8
    fingerprint    = '483448e46c2fe187'     ! optional ('' = none)
    species        = 'A', 'B', 'A', 'B', 'A', 'B', 'A', 'B'   ! optional
/

&REGION
    name  = 'top'
    atoms = 8
/

&REGION
    name  = 'interior'
    atoms = 2, 3, 4, 5, 6, 7
/
```

Rules (all violations are fatal and name the region and offending value):

- `&REGION_FILE` comes first, once; `format_version` must be 1; `natoms`
  must equal the structure `nat`.
- `species`, if present, must equal the structure labels in order.
- `fingerprint`, if non-empty, must equal the structure fingerprint.
- Any number of `&REGION` blocks. `name`: 1–64 characters from
  `[A-Za-z0-9_.]`, unique (case-sensitive), not `total` (reserved for the
  all-orbital group added by `spectral`). Names containing `+` or `-` are
  rejected, which keeps them addressable in `specmap` map expressions.
- `atoms`: non-empty, indices in `1..natoms`, no duplicates inside a region,
  at most 100000 entries per region.
- Regions may overlap (reported by `validate`). Overlapping regions are
  allowed but are double-counted if summed in `specmap`.
- Species labels longer than 8 characters are truncated by the namelist read
  (labels written by `atomselect` always fit).

### How a region is added to an existing file

`select`, `layers --regions` and `surface` all use the same logic:

1. If the regions file does **not** exist, a new one is started with a
   header taken from the structure: `format_version = 1`, `natoms`, the
   structure fingerprint and the species labels.
2. If it **exists**, it is read and validated against the structure
   (including the fingerprint; a file without a fingerprint is rejected
   here), the new region(s) are appended (the name must be unique), and the
   whole file is rewritten in canonical layout with the previous regions
   preserved in their order.

The file can also be written by hand; run `atomselect validate` on it.

## Fingerprint

The fingerprint ties a regions file to one structure (atom labels and
positions; the lattice vectors are **not** included). Canonical text, each
line terminated by a newline:

```text
LFSTRUCT1
<nat>
<label_1> <ix_1> <iy_1> <iz_1>
...
```

with `label = trim(adjustl(label))` and `ix = nint(x·10⁴)` written as a
plain integer (positions rounded to 10⁻⁴ Å). The fingerprint is the FNV-1a
64-bit hash of these bytes (offset basis `0xcbf29ce484222325`, prime
`0x100000001b3`), printed as 16 lowercase hex characters. Moving any atom by
more than ~10⁻⁴ Å, relabelling or reordering atoms changes it.

## Tests

- `fpm test test_regions` — `regions.nml` reading/writing/validation and
  numeric index lists.
- `fpm test test_selection_geom` — projection on arbitrary normals,
  range/box selection, layer detection, reference layers, surface layers,
  periodic-split warning.
- `fpm test test_structure_io` — structure reader and fingerprint (including
  the FNV-1a test vectors).
- `bash test/test_atomselect_cli.sh` — end-to-end CLI checks of every
  subcommand on `test/data/slab8.nml` (inspect/xyz/table, fingerprint
  stability, select by list/normal/box, overlaps and partition, layers with
  `--tol`/`--refs`/tilted normals, split warning, surface with `+z`, `-z` and
  tilted normals). Run from the repository root; it calls `fpm build`.

## Worked example

See [`examples/surface_spectral`](../examples/surface_spectral/README.md):

```bash
atomselect inspect  model.nml
atomselect layers   model.nml --normal 0 0 1
atomselect surface  model.nml --regions regions.nml --normal 0 0 1 --top 1 --bottom 1
atomselect select   model.nml --regions regions.nml --name interior --atoms "2-7"
atomselect validate model.nml regions.nml               # reports "partition   : yes"
```
