plot '../postProcessing/patchPanelQr/0/faceSource.dat' u 1:($3/1000) title "panel"
replot '../postProcessing/patchWallQr/0/faceSource.dat' u 1:($3/1000) title "wall"
replot '../postProcessing/patchSideQr/0/faceSource.dat' u 1:($3/1000) title "side"
replot '../postProcessing/patchOutletQr/0/faceSource.dat' u 1:($3/1000) title "outlet"
replot '../postProcessing/patchBurnerQr/0/faceSource.dat' u 1:($3/1000) title "burner"
replot '../postProcessing/patchGroundQr/0/faceSource.dat' u 1:($3/1000) title "ground"
replot "< paste ../postProcessing/patchPanelQr/0/faceSource.dat ../postProcessing/patchWallQr/0/faceSource.dat ../postProcessing/patchSideQr/0/faceSource.dat ../postProcessing/patchOutletQr/0/faceSource.dat ../postProcessing/patchBurnerQr/0/faceSource.dat ../postProcessing/patchGroundQr/0/faceSource.dat" u 1:(($3+$6+$9+$12+$15+$18)/1000) title "sum of above"
replot '../postProcessing/HRR/0/cellSource.dat' u 1:($3/1000) title "gas phase HRR"
replot '../postProcessing/HRR/0/cellSource.dat' u 1:($3/1000*0.6) title "HRR * RadFrac"
set xlabel "time [s]"
set ylabel "integrated heat flux [kW]"
set xrange [0:*]
set title "Integraed Net Radiative Heat Flux"
replot

