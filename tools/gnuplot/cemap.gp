# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Constant-energy map I(kx, ky) at a fixed energy (from a k-grid spectral file).
# Input: a specmap .cemap.dat file and the matching cegrid include.
#
# Usage:
#   gnuplot -e "datafile='slab.top.E+0.000.cemap.dat'; grid='slab.cegrid.gp'; out='slab_top_EF'" cemap.gp
# Optional:
#   title='...'
#   diff=1         divergent blue-white-red palette centred at 0 (for maps such as top-bottom)
#   cbmax=<value>  fixed colour limit ([0:cbmax], or [-cbmax:cbmax] with diff=1)
# The columns plotted come from the include (COLS = "1:2:6" kx,ky when the grid
# plane is perpendicular to z, otherwise "4:5:6" in-plane k_u,k_v).
# Writes <out>.png (pngcairo) and <out>.pdf (pdfcairo).

if (!exists("datafile")) { print "cemap.gp: set datafile='...cemap.dat'"; exit status 1 }
if (!exists("grid"))     { print "cemap.gp: set grid='...cegrid.gp'"; exit status 1 }
if (!exists("out"))      out = "cemap"
if (!exists("title"))    title = datafile
if (!exists("diff"))     diff = 0

# stats before loading the include: with a single column, stats filters the
# values through the current xrange, which the include sets to the k range.
stats datafile u 6 nooutput
load grid

set view map
set pm3d map
set pm3d corners2color c1
unset key
set border lw 1
set tics out nomirror
set title title noenhanced

if (diff) {
    m = (abs(STATS_min) > abs(STATS_max)) ? abs(STATS_min) : abs(STATS_max)
    if (m <= 0) { m = 1.0 }
    if (exists("cbmax")) { m = cbmax }
    set palette defined (-1 "#08306b", -0.5 "#4292c6", 0 "white", 0.5 "#ef3b2c", 1 "#67000d")
    set cbrange [-m:m]
    set cblabel "intensity difference (states/eV)"
} else {
    set palette defined (0 "black", 1 "#3b0f70", 2 "#8c2981", 3 "#de4968", 4 "#fe9f6d", 5 "#fcfdbf")
    set cbrange [0:*]
    if (exists("cbmax")) { set cbrange [0:cbmax] }
    set cblabel "intensity (states/eV)"
}

do for [fmt in "png pdf"] {
    if (fmt eq "png") { set terminal pngcairo size 900,800 enhanced font "Helvetica,14"; set output out.".png" }
    if (fmt eq "pdf") { set terminal pdfcairo size 4.5in,4in enhanced font "Helvetica,10"; set output out.".pdf" }
    splot datafile u @COLS w pm3d notitle
    set output
}
