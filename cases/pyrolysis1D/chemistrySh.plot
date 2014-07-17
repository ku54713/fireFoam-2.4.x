plot   './probes/panelRegion/0/chemistrySh' u 1:(-$2) title "1 cell"
replot './probes/panelRegion/0/chemistrySh' u 1:(-$3) title "2 cell"
replot './probes/panelRegion/0/chemistrySh' u 1:(-$4) title "3 cell"
replot './probes/panelRegion/0/chemistrySh' u 1:(-$5) title "4 cell"
replot './probes/panelRegion/0/chemistrySh' u 1:(-$6) title "5 cell"
replot './probes/panelRegion/0/chemistrySh' u 1:(-$7) title "6 cell"
replot './probes/panelRegion/0/chemistrySh' u 1:(-$8) title "7 cell"
replot './probes/panelRegion/0/chemistrySh' u 1:(-$9) title "8 cell"
set xlabel "Time [s]"
set ylabel "Heat absorption rate"
set xrange [0:*]
set title "Heat Absorption Rate"
set key right top
replot
set terminal postscript eps enhanced color
set output '| epstopdf --filter --outfile=plot.chemistrySh.pdf'
replot 
set output 
set terminal x11
