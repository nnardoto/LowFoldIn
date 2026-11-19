#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Testes de linha de comando do surfgreen (portão 2):
#  - execução com OpenMP (perfil release) e .green.bin idêntico bit a bit
#    para 1, 2 e 4 threads;
#  - --write-structure, --info e --slice (caminho k e grade k);
#  - códigos de saída em erro.
# Uso: bash test/test_surfgreen_cli.sh   (a partir da raiz do repositório)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
PROFILE="${PROFILE:-release}"
export OPENBLAS_NUM_THREADS=1
fpm build --profile "$PROFILE" >/dev/null
binpath() { echo "$ROOT/$(fpm run --profile "$PROFILE" "$@" --runner echo 2>/dev/null | tail -n 1)"; }
GEN="$(binpath --example surfgreen_model)"
SG="$(binpath surfgreen)"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
"$GEN" >/dev/null
pass=0; fail=0
ok()  { pass=$((pass + 1)); }
bad() { fail=$((fail + 1)); echo "  FAIL: $*"; }
check() { if eval "$1"; then ok; else bad "$2"; fi; }

cat > s.nml <<'NML'
&GREEN
    mode = 'film', eta_ev = 0.02, emin = -2.0, emax = 2.0, ne = 31, output_prefix = 'run'
/
&SUBSTRATE
    h_file = 'sub_H.hsio', s_file = 'sub_S.hsio', structure_file = 'sub_struct.nml', central_layers = 1
/
&FILM
    h_file = 'film_H.hsio', s_file = 'film_S.hsio', structure_file = 'film_struct.nml', nlayers = 3
/
&INTERFACE
    model = 'scaled', lambda_h = 0.6, lambda_s = 0.6
/
&COHERENT
    nqz = 2, qz_min = 0.0, qz_max = 1.0, escape_length_ang = 5.0
/
&KPATH
    nvertex = 2
    vertex(:,1) = 0.0, 0.0, 0.0
    vertex(:,2) = 0.5, 0.0, 0.0
    label = 'G', 'X'
    npoints = 12
/
NML

for t in 1 2 4; do
    OMP_NUM_THREADS=$t "$SG" s.nml >/dev/null
    mv run.green.bin "run_t$t.green.bin"
done
check "cmp -s run_t1.green.bin run_t2.green.bin" "green.bin differs between 1 and 2 threads"
check "cmp -s run_t1.green.bin run_t4.green.bin" "green.bin differs between 1 and 4 threads"
check "[[ -s run_central.nml ]]" "central structure not written"

"$SG" --info run_t1.green.bin > info.txt
check "grep -q 'mode = film' info.txt" "--info lacks mode"
check "grep -q 'L004_film' info.txt" "--info lacks default layer regions"

"$SG" --slice run_t1.green.bin --map L002_film --out sl >/dev/null
check "[[ \$(grep -vc '^#' sl.map.dat) -eq \$((13 * 31 + 13)) ]]" "path slice has wrong number of lines"
check "grep -q 'set xtics' sl.kticks.gp" "kticks include lacks xtics"
"$SG" --slice run_t1.green.bin --map coh --iq 2 --out co >/dev/null
check "[[ -s co.map.dat ]]" "coherent slice not written"

rm -f run_central.nml
"$SG" --write-structure s.nml >/dev/null
check "[[ -s run_central.nml && ! -e run.green.bin ]]" "--write-structure should only write the structure"

sed -e "s/ne = 31/ne = 1/" -e "s/emax = 2.0, //" s.nml > g.nml
python_free_grid=$'&KGRID\n    origin = -0.2, -0.2, 0.0\n    span1 = 0.4, 0.0, 0.0\n    span2 = 0.0, 0.4, 0.0\n    n1 = 5, n2 = 4\n/'
awk '/&KPATH/{exit} {print}' g.nml > g2.nml
printf '%s\n' "$python_free_grid" >> g2.nml
"$SG" g2.nml >/dev/null
"$SG" --slice run.green.bin --out ce >/dev/null
check "[[ \$(grep -vc '^#' ce.cemap.dat) -eq \$((5 * 4 + 4)) ]]" "grid slice has wrong number of lines"
check "grep -q 'COLS = \"1:2:6\"' ce.cegrid.gp" "cegrid include lacks COLS"

if "$SG" missing.nml 2>/dev/null; then bad "missing input should fail"; else ok; fi
sed "s/model = 'scaled'/model = 'extracted'/" s.nml > e.nml
if "$SG" e.nml 2>err.txt; then bad "'extracted' should fail"; else ok; fi
check "grep -q 'phase 4' err.txt" "'extracted' error message"

echo "[test_surfgreen_cli] $pass passed, $fail failed"
[[ $fail -eq 0 ]]
