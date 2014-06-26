plot '../postProcessing/patchPanelQr/0/faceSource.dat' u 1:($3/1000) title "Net Radiative"
replot '../postProcessing/patchPanelConvectiveHeatFlux/0/faceSource.dat' u 1:($3/1000) title "Convective"
set xlabel "time [s]"
set ylabel "integrated heat flux [kW]"
set xrange [0:*]
set title "Integraed Heat Flux"
replot

