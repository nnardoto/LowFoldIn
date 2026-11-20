#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Slater-Koster contact on SYNTHETIC materials (NOT DFT):
#   gen_materials -> skfit (A) -> skfit (B) -> skfit (combine) -> surfgreen &SCAN
#   -> surfgreen --slice -> gnuplot
#
# Usage: bash examples/sk_contact/run_example.sh
#        NO_PLOT=1 bash examples/sk_contact/run_example.sh
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
PROFILE="${PROFILE:-release}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

cd "$ROOT"
fpm build --profile "$PROFILE"
binpath() { echo "$ROOT/$(fpm run --profile "$PROFILE" "$@" --runner echo 2>/dev/null | tail -n 1)"; }
GEN="$(binpath --example sk_materials)"
SKFIT="$(binpath skfit)"
SURFGREEN="$(binpath surfgreen)"

cd "$HERE"
step() { printf '\n==== %s ====\n' "$*"; }

step "0. clean"
rm -f ./*.hsio A.nml B.nml ./*.sk ./*_central.nml ./*.green.bin ./*.scan.txt ./*.map.dat ./*.kticks.gp ./*.png ./*.pdf

step "1. synthetic materials (NOT DFT)"
"$GEN"

step "2. fit SK to the coupling between principal layers of each bulk"
"$SKFIT" fit_A.nml
"$SKFIT" fit_B.nml

step "3. contact A-B by the combination rule"
"$SKFIT" combine.nml
cat contact.sk

step "4. film of B on semi-infinite A, scan of lambda_h"
"$SURFGREEN" film_scan.nml
cat b_on_a.scan.txt
for f in b_on_a_s00?.green.bin; do
    "$SURFGREEN" --slice "$f" --map total --out "${f%.green.bin}.total"
done

step "5. gnuplot"
GP="$ROOT/tools/gnuplot"
if [[ "${NO_PLOT:-0}" == "1" ]]; then
    echo "NO_PLOT=1: skipping gnuplot step"
elif ! command -v gnuplot >/dev/null 2>&1; then
    echo "gnuplot not found: skipping plot step"
else
    for f in b_on_a_s00?.total.map.dat; do
        b="${f%.map.dat}"
        gnuplot -e "datafile='$f'; ticks='$b.kticks.gp'; out='$b'" "$GP/map_positive.gp"
    done
    ls -1 ./*.png
fi
step "done"
