#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Full surface-spectral pipeline on a SYNTHETIC model slab (NOT DFT, NOT ARPES):
#   gen_model -> atomselect -> spectral -> specmap (x2) -> gnuplot
#
# Usage (from anywhere; the script cd's into its own directory):
#   bash examples/surface_spectral/run_example.sh
#   NO_PLOT=1 bash examples/surface_spectral/run_example.sh   # skip gnuplot
# The gnuplot step is also skipped when gnuplot is not installed.
# Binaries are built with fpm (default profile) from the repository root.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

# ---- build and locate binaries ----------------------------------------------
cd "$ROOT"
fpm build
binpath() { echo "$ROOT/$(fpm run "$@" --runner echo 2>/dev/null | tail -n 1)"; }
GEN="$(binpath --example surface_model)"
ATOMSELECT="$(binpath atomselect)"
SPECTRAL="$(binpath spectral)"
SPECMAP="$(binpath specmap)"
for b in "$GEN" "$ATOMSELECT" "$SPECTRAL" "$SPECMAP"; do
    [[ -x "$b" ]] || { echo "binary not found: $b" >&2; exit 1; }
done

cd "$HERE"
step() { printf '\n==== %s ====\n' "$*"; }

# ---- 0. clean previous outputs (regions.nml is regenerated from scratch) ----
step "0. clean"
rm -f model_H.hsio model_S.hsio model.nml regions.nml \
      ./*.spectral.bin ./*.spectral.txt ./*.map.dat ./*.kticks.gp ./*.png ./*.pdf

# ---- 1. synthetic model -----------------------------------------------------
step "1. generate synthetic model (NOT DFT)"
"$GEN"

# ---- 2. regions ---------------------------------------------------------------
step "2a. atomselect inspect"
"$ATOMSELECT" inspect model.nml
step "2b. atomselect layers"
"$ATOMSELECT" layers model.nml --normal 0 0 1
step "2c. atomselect surface (top/bottom)"
"$ATOMSELECT" surface model.nml --regions regions.nml --normal 0 0 1 --top 1 --bottom 1
step "2d. atomselect select (interior)"
"$ATOMSELECT" select model.nml --regions regions.nml --name interior --atoms "2-7"

# ---- 3. validate ----------------------------------------------------------------
step "3. atomselect validate"
"$ATOMSELECT" validate model.nml regions.nml

# ---- 4. spectral ----------------------------------------------------------------
step "4. spectral (Mulliken, X-G-M)"
"$SPECTRAL" spectral.nml

# ---- 5./6. specmap ----------------------------------------------------------------
step "5. specmap (Lorentzian, eta = 0.05 eV)"
"$SPECMAP" specmap.nml
step "6. specmap (Gaussian, sigma = 0.10 eV; same spectral file, no diagonalisation)"
"$SPECMAP" specmap_gauss.nml

# ---- 7. plots -------------------------------------------------------------------
step "7. gnuplot"
if [[ "${NO_PLOT:-0}" == "1" ]]; then
    echo "NO_PLOT=1: skipping gnuplot step"
elif ! command -v gnuplot >/dev/null 2>&1; then
    echo "gnuplot not found: skipping plot step"
else
    gnuplot plot_spectral.gp
    ls -1 model_*.png model_*.pdf
fi

step "done"
