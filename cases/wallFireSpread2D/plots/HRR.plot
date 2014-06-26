plot '../postProcessing/HRR/0/cellSource.dat' u 1:($3/1000) title "gas phase"
replot '../postProcessing/patchPanel/0/faceSource.dat' u 1:($4/1000) title "panel"
replot '../postProcessing/patchBurner/0/faceSource.dat' u 1:($4/1000) title "burner"
replot '../postProcessing/patchOutlet/0/faceSource.dat' u 1:(-$4/1000) title "outlet"
replot "< paste ../postProcessing/HRR/0/cellSource.dat ../postProcessing/patchOutlet/0/faceSource.dat" u 1:(($3-$7)/1000) title "gas phase + outlet"
replot "< paste ../postProcessing/patchPanel/0/faceSource.dat ../postProcessing/patchBurner/0/faceSource.dat" u 1:(($4+$8)/1000) title "burner + panel"
set xlabel "time (s)"
set ylabel "HRR [kW]"
set xrange [0:*]
set title "Heat releast rate"
set key left top
replot

set terminal postscript eps enhanced color
set output '| epstopdf --filter --outfile=plot.HRR.pdf'
replot 
set output 
set terminal x11
