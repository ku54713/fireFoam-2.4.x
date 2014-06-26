plot '../postProcessing/patchPanel/0/faceSource.dat' u 1:(-$3) title "phi (gas region)"
replot '../postProcessing/patchPanelSolidRegion/0/faceSource.dat' u 1:3 title "phiGas (solid region)"
set xlabel "time [s]"
set ylabel "phi [kg/s]"
set xrange [0:*]
set title "fuel mass flow rate"
set key left top
replot

set terminal postscript eps enhanced color
set output '| epstopdf --filter --outfile=plot.phiPanel.pdf'
replot 
set output 
set terminal x11
