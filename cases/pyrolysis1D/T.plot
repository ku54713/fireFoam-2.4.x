plot   './probes/panelRegion/0/T' u 1:2 title "1 cell"
replot './probes/panelRegion/0/T' u 1:3 title "2 cell"
replot './probes/panelRegion/0/T' u 1:4 title "3 cell"
replot './probes/panelRegion/0/T' u 1:5 title "4 cell"
replot './probes/panelRegion/0/T' u 1:6 title "5 cell"
replot './probes/panelRegion/0/T' u 1:7 title "6 cell"
replot './probes/panelRegion/0/T' u 1:8 title "7 cell"
replot './probes/panelRegion/0/T' u 1:9 title "8 cell"
replot './patchProbes/panelRegion/0/T' u 1:2 title "surface"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set xrange [0:*]
set title "Solid Phase Temperature"
set key right bottom
replot
set terminal postscript eps enhanced color
set output '| epstopdf --filter --outfile=plot.T.pdf'
replot 
set output 
set terminal x11
