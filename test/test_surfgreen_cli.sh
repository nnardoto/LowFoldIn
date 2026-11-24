#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Testes de linha de comando do surfgreen (portão 2):
#  - execução com OpenMP (perfil release) e .green.bin idêntico bit a bit
#    para 1, 2 e 4 threads;
#  - --write-structure, --info e --slice (caminho k e grade k);
#  - &ARPES, --kperp-slice e reprodutibilidade do eixo físico k_perp;
#  - --photon-slice, conversão He-I, referências de energia e threads;
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

sed -e "s/&COHERENT/\&ARPES/" -e "s/nqz = 2/nkperp = 2/" \
    -e "s/qz_min = 0.0/qz_min = 2.0/" -e "s/qz_max = 1.0/qz_max = 3.0/" \
    -e "s/qz_min/kperp_min/" -e "s/qz_max/kperp_max/" \
    -e "s/output_prefix = 'run'/output_prefix = 'arpes'/" s.nml > a.nml
for t in 1 2 4; do
    OMP_NUM_THREADS=$t "$SG" a.nml >/dev/null
    mv arpes.green.bin "arpes_t$t.green.bin"
done
check "cmp -s arpes_t1.green.bin arpes_t2.green.bin" "ARPES green.bin differs between 1 and 2 threads"
check "cmp -s arpes_t1.green.bin arpes_t4.green.bin" "ARPES green.bin differs between 1 and 4 threads"
"$SG" --info arpes_t1.green.bin > arpes_info.txt
check "grep -q 'detection_axis = k_perp' arpes_info.txt" "ARPES --info lacks physical axis metadata"
check "grep -q 'k_perp range' arpes_info.txt" "ARPES --info lacks k_perp range"
"$SG" --kperp-slice arpes_t1.green.bin --ik 2 --out ak >/dev/null
check "[[ \$(grep -vc '^#' ak.kperp.dat) -eq \$((2 * 31 + 2)) ]]" "k_perp slice has wrong number of lines"
check "grep -q '# columns: k_perp_invA energy_eV intensity' ak.kperp.dat" "k_perp slice lacks units"
"$SG" --slice arpes_t1.green.bin --map coh --iq 2 --out ac >/dev/null
check "grep -q '# k_perp = ' ac.map.dat" "ARPES k-path slice lacks selected k_perp"
if "$SG" --kperp-slice run_t1.green.bin >/dev/null 2>&1; then
    bad "legacy q_z file should not have an ARPES k_perp axis"
else
    ok
fi

"$SG" --photon-slice arpes_t1.green.bin --photon-energy 21.2 --work-function 4.5 \
    --inner-potential 10.0 --energy-reference E-E_F --ik 1 --out ph >/dev/null
check "[[ \$(grep -vc '^#' ph.photon.dat) -eq 31 ]]" "photon slice has wrong number of energies"
check "grep -q '# columns: initial_energy_eV kinetic_energy_eV k_perp_invA intensity' ph.photon.dat" \
    "photon slice lacks physical columns"
check "awk '!/^#/{e=\$1; ek=\$2; k=\$3; exit} END{exit !(e+2.0 < 1e-8 && e+2.0 > -1e-8 && ek-14.7 < 1e-8 && ek-14.7 > -1e-8 && k-2.546167721 < 1e-8 && k-2.546167721 > -1e-8)}' ph.photon.dat" \
    "He-I conversion at E=-2 eV disagrees with analytic value"
"$SG" --kperp-slice arpes_t1.green.bin --ik 1 --out akg >/dev/null
check "awk 'FNR==NR{if (\$0 !~ /^#/ && NF==3){n++; if(n==1)i1=\$3; if(n==32)i2=\$3} next} !/^#/{x=0.54616772093; expected=(1-x)*i1+x*i2; exit !(\$4-expected < 1e-8 && \$4-expected > -1e-8)}' akg.kperp.dat ph.photon.dat" \
    "photon slice intensity is not the independent linear k_perp interpolation"
"$SG" --photon-slice arpes_t4.green.bin --photon-energy 21.2 --work-function 4.5 \
    --inner-potential 10.0 --energy-reference E-E_F --ik 1 --out ph4 >/dev/null
check "diff -q <(grep -v '^#' ph.photon.dat) <(grep -v '^#' ph4.photon.dat) >/dev/null" \
    "photon slice differs between 1 and 4 threads"
"$SG" --photon-slice arpes_t1.green.bin --photon-energy 21.2 --work-function 4.5 \
    --inner-potential 10.0 --energy-reference absolute --fermi-energy -4.25 --ik 1 --out pha >/dev/null
check "awk '!/^#/{exit !(\$2-18.95 < 1e-8 && \$2-18.95 > -1e-8)}' pha.photon.dat" \
    "absolute energy reference does not use E-E_F"
if "$SG" --photon-slice arpes_t1.green.bin --photon-energy 21.2 --work-function 4.5 \
    --inner-potential 10.0 --energy-reference absolute --ik 1 >/dev/null 2>&1; then
    bad "absolute film energies should require explicit E_F"
else
    ok
fi
if "$SG" --photon-slice arpes_t1.green.bin --photon-energy 21.2 --work-function 4.5 \
    --inner-potential 10.0 --energy-reference E-E_F --ik 2 >/dev/null 2>&1; then
    bad "off-normal k-point should not use the normal-emission photon conversion"
else
    ok
fi
if "$SG" --photon-slice arpes_t1.green.bin --photon-energy 21.2 --work-function 4.5 \
    --energy-reference E-E_F >/dev/null 2>&1; then
    bad "photon slice should require explicit V0"
else
    ok
fi

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

# modo 'extracted' (fase 4): o substrato faz o papel da supercélula e do lead
cat > x.nml <<'NML'
&GREEN
    mode = 'extracted', eta_ev = 0.02, emin = -2.0, emax = 2.0, ne = 31, output_prefix = 'xrun'
/
&EXTRACTED
    h_file = 'sub_H.hsio', s_file = 'sub_S.hsio', structure_file = 'sub_struct.nml'
    embed_axis = 3, central_columns = 0, 2, left = 'lead', right = 'lead'
/
&LEAD
    h_file = 'sub_H.hsio', s_file = 'sub_S.hsio', structure_file = 'sub_struct.nml'
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
    OMP_NUM_THREADS=$t "$SG" x.nml >/dev/null
    mv xrun.green.bin "xrun_t$t.green.bin"
done
check "cmp -s xrun_t1.green.bin xrun_t2.green.bin" "extracted green.bin differs between 1 and 2 threads"
check "cmp -s xrun_t1.green.bin xrun_t4.green.bin" "extracted green.bin differs between 1 and 4 threads"
"$SG" --info xrun_t1.green.bin > xinfo.txt
check "grep -q 'mode = extracted' xinfo.txt" "extracted --info lacks mode"
check "grep -q 'col1' xinfo.txt" "extracted --info lacks default column regions"
"$SG" --slice xrun_t1.green.bin --map col1 --out xs >/dev/null
check "[[ \$(grep -vc '^#' xs.map.dat) -eq \$((13 * 31 + 13)) ]]" "extracted path slice has wrong number of lines"
rm -f xrun_central.nml
"$SG" --write-structure x.nml >/dev/null
check "[[ -s xrun_central.nml && ! -e xrun.green.bin ]]" "extracted --write-structure should only write the structure"
sed "s/right = 'lead'/right = 'vacuum'/" x.nml > xe.nml
if "$SG" xe.nml 2>xerr.txt; then bad "extracted with atoms beyond a vacuum side should fail"; else ok; fi
check "grep -q 'neither central nor' xerr.txt" "extracted error message for an unmatched coupling"

echo "[test_surfgreen_cli] $pass passed, $fail failed"
[[ $fail -eq 0 ]]
