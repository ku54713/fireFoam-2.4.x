plot '../postProcessing/patchProbes/0/Qin' u 1:($2/1000) title "0.225"
replot '../postProcessing/patchProbes/0/Qin' u 1:($3/1000) title "0.325"
replot '../postProcessing/patchProbes/0/Qin' u 1:($4/1000) title "0.425"
replot '../postProcessing/patchProbes/0/Qin' u 1:($5/1000) title "0.525"
set xlabel "time [s]"
set ylabel "Qr [kW/m^2]"
set xrange [0:*]
set title "Incident Radiative Heat Flux"
replot

