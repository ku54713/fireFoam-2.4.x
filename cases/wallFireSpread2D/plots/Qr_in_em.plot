plot '../postProcessing/patchProbes/0/Qr' u 1:($2/1000) title "Qr, 0.225"
#replot '../postProcessing/patchProbes/0/Qr' u 1:($3/1000) title "Qr, 0.325"
#replot '../postProcessing/patchProbes/0/Qr' u 1:($4/1000) title "Qr, 0.425"
#replot '../postProcessing/patchProbes/0/Qr' u 1:($5/1000) title "Qr, 0.525"
replot '../postProcessing/patchProbes/0/Qin' u 1:($2/1000) title "Qin, 0.225"
#replot '../postProcessing/patchProbes/0/Qin' u 1:($3/1000) title "Qin, 0.325"
#replot '../postProcessing/patchProbes/0/Qin' u 1:($4/1000) title "Qin, 0.425"
#replot '../postProcessing/patchProbes/0/Qin' u 1:($5/1000) title "Qin, 0.525"
replot "< paste ../postProcessing/patchProbes/0/Qin ../postProcessing/patchProbes/0/Qem" u 1:(($2+$7)/1000) title "Qin+Qem, 0.225"
#replot "< paste ../postProcessing/patchProbes/0/Qin ../postProcessing/patchProbes/0/Qem" u 1:(($3+$8)/1000) title "Qin+Qem, 0.325"
#replot "< paste ../postProcessing/patchProbes/0/Qin ../postProcessing/patchProbes/0/Qem" u 1:(($4+$9)/1000) title "Qin+Qem, 0.425"
#replot "< paste ../postProcessing/patchProbes/0/Qin ../postProcessing/patchProbes/0/Qem" u 1:(($5+$10)/1000) title "Qin+Qem, 0.525"
set xlabel "time [s]"
set ylabel "Qr [kW/m^2]"
set xrange [0:*]
set title "Incident Radiative Heat Flux"
replot

