# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Side-by-side comparison (multiplot 1x2) of two spectral maps sharing the
# same k-path, with the sequential palette and a common colour range.
#
# Usage:
#   gnuplot -e "datafile1='slab.top.map.dat'; datafile2='slab.bottom.map.dat'; ticks='slab.kticks.gp'; out='slab_cmp'" map_compare.gp
# Optional: title1='...', title2='...', cbmax=<value>.
# Writes <out>.png (pngcairo) and <out>.pdf (pdfcairo).

if (!exists("datafile1") || !exists("datafile2")) { print "map_compare.gp: set datafile1 and datafile2"; exit status 1 }
if (!exists("ticks"))  { print "map_compare.gp: set ticks='...kticks.gp'"; exit status 1 }
if (!exists("out"))    out = "map_compare"
if (!exists("title1")) title1 = datafile1
if (!exists("title2")) title2 = datafile2

VLINE_COLOR = "white"
load ticks

# Common colour range: [0 : max of both maps]
stats datafile1 u 3 nooutput
m1 = STATS_max
stats datafile2 u 3 nooutput
m2 = STATS_max
m = (m1 > m2) ? m1 : m2
if (m <= 0) m = 1.0
if (exists("cbmax")) m = cbmax
stats datafile1 u 2 nooutput
set yrange [STATS_min:STATS_max]

set view map
set pm3d map
set pm3d corners2color c1
unset key
set tics out nomirror
set xtics scale 0
set xlabel "k"
set palette defined (0 "black", 1 "#3b0f70", 2 "#8c2981", 3 "#de4968", 4 "#fe9f6d", 5 "#fcfdbf")
set cbrange [0:m]
set cblabel "intensity (states/eV)"

if (EF == EF) set arrow from graph 0, first EF to graph 1, first EF nohead dt 2 lw 1.5 lc rgb "white" front

do for [fmt in "png pdf"] {
    if (fmt eq "png") { set terminal pngcairo size 1600,700 enhanced font "Helvetica,14"; set output out.".png" }
    if (fmt eq "pdf") { set terminal pdfcairo size 9in,4in enhanced font "Helvetica,10"; set output out.".pdf" }
    set multiplot layout 1,2
    set ylabel ELABEL offset -1,0
    set title title1 noenhanced
    unset colorbox
    splot datafile1 u 1:2:3 w pm3d notitle
    set ylabel ""
    set title title2 noenhanced
    set colorbox
    splot datafile2 u 1:2:3 w pm3d notitle
    unset multiplot
    set output
}
