# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Example-specific wrapper around tools/gnuplot/*.gp for the SYNTHETIC model
# slab (NOT DFT, NOT ARPES). Run from examples/surface_spectral after specmap:
#
#   gnuplot plot_spectral.gp
#
# Writes PNG + PDF:
#   model_total      total A(k,E)                          (map_positive.gp)
#   model_top        top-layer projected A(k,E)            (map_positive.gp)
#   model_top_minus_bottom  top - bottom, divergent palette (map_difference.gp)
#   model_top_vs_bottom     top | bottom, common colour scale (map_compare.gp)
# Each block starts with "reset session" so that the generic scripts never
# inherit variables or settings from the previous plot.

reset session
datafile = 'model.total.map.dat'; ticks = 'model.kticks.gp'; out = 'model_total'
title = 'total (synthetic model, Lorentzian HWHM 0.05 eV)'
load '../../tools/gnuplot/map_positive.gp'

reset session
datafile = 'model.top.map.dat'; ticks = 'model.kticks.gp'; out = 'model_top'
title = 'top layer (atom 8) (synthetic model)'
load '../../tools/gnuplot/map_positive.gp'

reset session
datafile = 'model.top_minus_bottom.map.dat'; ticks = 'model.kticks.gp'; out = 'model_top_minus_bottom'
title = 'top - bottom (red: top, blue: bottom) (synthetic model)'
load '../../tools/gnuplot/map_difference.gp'

reset session
datafile1 = 'model.top.map.dat'; datafile2 = 'model.bottom.map.dat'
ticks = 'model.kticks.gp'; out = 'model_top_vs_bottom'
title1 = 'top layer (atom 8)'; title2 = 'bottom layer (atom 1)'
load '../../tools/gnuplot/map_compare.gp'
