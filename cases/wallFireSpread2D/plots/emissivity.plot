plot '../postProcessing/patchProbes/panelRegion/0/emissivity' u 1:2 title "0.225"
replot '../postProcessing/patchProbes/panelRegion/0/emissivity' u 1:3 title "0.325"
replot '../postProcessing/patchProbes/panelRegion/0/emissivity' u 1:4 title "0.425"
replot '../postProcessing/patchProbes/panelRegion/0/emissivity' u 1:5 title "0.525"
set xlabel "time (s)"
set ylabel "emissivity"
set xrange [0:*]
set title "emissivity"
replot

