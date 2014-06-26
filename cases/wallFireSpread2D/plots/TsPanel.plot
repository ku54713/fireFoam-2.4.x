plot '../postProcessing/patchProbes/panelRegion/0/T' u 1:2 title "0.225"
replot '../postProcessing/patchProbes/panelRegion/0/T' u 1:3 title "0.325"
replot '../postProcessing/patchProbes/panelRegion/0/T' u 1:4 title "0.425"
replot '../postProcessing/patchProbes/panelRegion/0/T' u 1:5 title "0.525"
set xlabel "time [s]"
set ylabel "T [K]"
set xrange [0:*]
set title "Panel Surface Temperature"
set key left top
replot

set terminal postscript eps enhanced color
set output '| epstopdf --filter --outfile=plot.TsPanel.pdf'
replot 
set output 
set terminal x11
