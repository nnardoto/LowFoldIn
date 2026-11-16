# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Spectral map A_P(k,E) with a sequential palette (non-negative intensities).
# Input: a specmap map file and the matching kticks include.
#
# Usage:
#   gnuplot -e "datafile='slab.top.map.dat'; ticks='slab.kticks.gp'; out='slab_top'" map_positive.gp
# Optional: title='...', cbmax=<value> (fixed upper colour limit).
# Writes <out>.png (pngcairo) and <out>.pdf (pdfcairo).

if (!exists("datafile")) { print "map_positive.gp: set datafile='...'"; exit status 1 }
if (!exists("ticks"))    { print "map_positive.gp: set ticks='...kticks.gp'"; exit status 1 }
if (!exists("out"))      out = "map_positive"
if (!exists("title"))    title = datafile

VLINE_COLOR = "white"
load ticks

set view map
set pm3d map
set pm3d corners2color c1
unset key
set border lw 1
set tics out nomirror
set xtics scale 0
set ylabel ELABEL offset -1,0
set xlabel "k"
set cblabel "intensity (states/eV)"
set title title noenhanced
set palette defined (0 "black", 1 "#3b0f70", 2 "#8c2981", 3 "#de4968", 4 "#fe9f6d", 5 "#fcfdbf")
set cbrange [0:*]
if (exists("cbmax")) set cbrange [0:cbmax]
stats datafile u 2 nooutput
set yrange [STATS_min:STATS_max]

# Fermi level (EF is NaN when unknown in the file's reference; NaN != NaN)
if (EF == EF) set arrow from graph 0, first EF to graph 1, first EF nohead dt 2 lw 1.5 lc rgb "white" front

do for [fmt in "png pdf"] {
    if (fmt eq "png") { set terminal pngcairo size 1000,800 enhanced font "Helvetica,14"; set output out.".png" }
    if (fmt eq "pdf") { set terminal pdfcairo size 5in,4in enhanced font "Helvetica,10"; set output out.".pdf" }
    splot datafile u 1:2:3 w pm3d notitle
    set output
}
