#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# surfgreen pipeline on SYNTHETIC models (NOT DFT):
#   gen_model -> surfgreen (surface) -> surfgreen --write-structure -> atomselect
#   -> surfgreen (film, cemap + coherent) -> surfgreen --slice -> gnuplot
#
# Usage: bash examples/surfgreen_model/run_example.sh
#        NO_PLOT=1 bash examples/surfgreen_model/run_example.sh
# PROFILE=release (default) builds with OpenMP; keep BLAS single-threaded.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
PROFILE="${PROFILE:-release}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

cd "$ROOT"
fpm build --profile "$PROFILE"
binpath() { echo "$ROOT/$(fpm run --profile "$PROFILE" "$@" --runner echo 2>/dev/null | tail -n 1)"; }
GEN="$(binpath --example surfgreen_model)"
SURFGREEN="$(binpath surfgreen)"
ATOMSELECT="$(binpath atomselect)"

cd "$HERE"
step() { printf '\n==== %s ====\n' "$*"; }

step "0. clean"
rm -f ./*.hsio sub_struct.nml film_struct.nml ./*_central.nml regions.nml ./*.green.bin ./*.map.dat ./*.kticks.gp \
      ./*.cemap.dat ./*.cegrid.gp ./*.png ./*.pdf

step "1. synthetic models (NOT DFT)"
"$GEN"

step "2. semi-infinite substrate surface"
"$SURFGREEN" surface.nml
"$SURFGREEN" --slice sub_surface.green.bin --map L003_sub --out sub_surface.L3

step "3. central structure of the film system and regions"
"$SURFGREEN" --write-structure film.nml
"$ATOMSELECT" surface film_on_sub_central.nml --regions regions.nml --normal 0 0 1 --top 1 --bottom 1
"$ATOMSELECT" select film_on_sub_central.nml --regions regions.nml --name film_top --atoms "5"
"$ATOMSELECT" validate film_on_sub_central.nml regions.nml

step "4. film on substrate (E x k)"
"$SURFGREEN" film.nml
"$SURFGREEN" --info film_on_sub.green.bin
"$SURFGREEN" --slice film_on_sub.green.bin --map top --out film_on_sub.top
"$SURFGREEN" --slice film_on_sub.green.bin --map total --out film_on_sub.total

step "5. constant-energy map and coherent q_z scan"
"$SURFGREEN" cemap.nml
"$SURFGREEN" --slice film_cemap.green.bin --map total --out film_cemap.total
"$SURFGREEN" --slice film_cemap.green.bin --map coh --iq 2 --out film_cemap.coh2

step "6. gnuplot"
GP="$ROOT/tools/gnuplot"
if [[ "${NO_PLOT:-0}" == "1" ]]; then
    echo "NO_PLOT=1: skipping gnuplot step"
elif ! command -v gnuplot >/dev/null 2>&1; then
    echo "gnuplot not found: skipping plot step"
else
    gnuplot -e "datafile='sub_surface.L3.map.dat'; ticks='sub_surface.L3.kticks.gp'; out='sub_surface_L3'" "$GP/map_positive.gp"
    gnuplot -e "datafile='film_on_sub.top.map.dat'; ticks='film_on_sub.top.kticks.gp'; out='film_on_sub_top'" "$GP/map_positive.gp"
    gnuplot -e "datafile='film_cemap.total.cemap.dat'; grid='film_cemap.total.cegrid.gp'; out='film_cemap_total'" "$GP/cemap.gp"
    ls -1 ./*.png
fi
step "done"
