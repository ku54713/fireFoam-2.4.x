plot '../postProcessing/patchProbes/0/T' u 1:2 title "gas phase (surface)"
replot '../postProcessing/patchProbes/panelRegion/0/T' u 1:2 title "solid phase (surface)"
replot '../postProcessing/probes/0/T' u 1:2 title "gas phase (internal)"
replot '../postProcessing/probes/panelRegion/0/T' u 1:2 title "solid phase (internal)"
set xlabel "time (s)"
set ylabel "Temperature [K]"
set xrange [0:*]
set title "Temperature"
replot

