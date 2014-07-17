plot './patchProbes/0/T' u 1:2 title "gas phase"
replot './patchProbes/panelRegion/0/T' u 1:2 title "solid phase"
set xlabel "time (s)"
set ylabel "Temperature [K]"
set xrange [0:*]
set title "Surface Temperature"
replot

