clear
plot './patchPanelSolid/0/faceSource.dat' u 1:3 title "solid phase patch"
replot './patchPanel/0/faceSource.dat' u 1:(-$3) title "gas phase patch"
replot './referenceResult/patchPanelSolid/0/faceSource.dat' u 1:3 title "reference"
set xlabel "Time [s]"
set ylabel "fuel mass loss rate [kg/s]"
set xrange [0:*]
set title "Fuel Mass Loss Rate"
set key right top
replot

set terminal postscript eps enhanced color
set output '| epstopdf --filter --outfile=plot.mlr.pdf'
replot 
set output 
set terminal x11
