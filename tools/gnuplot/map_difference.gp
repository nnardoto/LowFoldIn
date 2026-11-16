# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Difference map (e.g. top-bottom) with a divergent blue-white-red palette
# centred at zero: cbrange = [-m:m], m = max |intensity| (from stats).
#
# Usage:
#   gnuplot -e "datafile='slab.top_minus_bottom.map.dat'; ticks='slab.kticks.gp'; out='slab_diff'" map_difference.gp
# Optional: title='...', cbmax=<m> (fixed symmetric colour limit).
# Writes <out>.png (pngcairo) and <out>.pdf (pdfcairo).

if (!exists("datafile")) { print "map_difference.gp: set datafile='...'"; exit status 1 }
if (!exists("ticks"))    { print "map_difference.gp: set ticks='...kticks.gp'"; exit status 1 }
if (!exists("out"))      out = "map_difference"
if (!exists("title"))    title = datafile

VLINE_COLOR = "gray40"
load ticks

stats datafile u 3 nooutput
m = (abs(STATS_min) > abs(STATS_max)) ? abs(STATS_min) : abs(STATS_max)
if (m <= 0) m = 1.0
if (exists("cbmax")) m = cbmax
stats datafile u 2 nooutput
set yrange [STATS_min:STATS_max]

set view map
set pm3d map
set pm3d corners2color c1
unset key
set tics out nomirror
set xtics scale 0
set ylabel ELABEL offset -1,0
set xlabel "k"
set cblabel "intensity difference (states/eV)"
set title title noenhanced
set palette defined (-1 "#08306b", -0.5 "#4292c6", 0 "white", 0.5 "#ef3b2c", 1 "#67000d")
set cbrange [-m:m]

if (EF == EF) set arrow from graph 0, first EF to graph 1, first EF nohead dt 2 lw 1.5 lc rgb "black" front

do for [fmt in "png pdf"] {
    if (fmt eq "png") { set terminal pngcairo size 1000,800 enhanced font "Helvetica,14"; set output out.".png" }
    if (fmt eq "pdf") { set terminal pdfcairo size 5in,4in enhanced font "Helvetica,10"; set output out.".pdf" }
    splot datafile u 1:2:3 w pm3d notitle
    set output
}
