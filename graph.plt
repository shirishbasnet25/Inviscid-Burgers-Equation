set multiplot layout 1,3

set grid 
set key off 

set size ratio 0.75
set title "{/:Bold Lax Method to Inviscid Burgers Equation}"
set xlabel ("x (m)")
set ylabel ("u (m/s)")
set xrange[0:1]
set yrange[-2:2]
plot "output.txt" using 1:2 w l lw 1.5 lc rgb 'green'
set title "{/:Bold Lax Wendroff Method to Inviscid Burgers Equation}"
plot "output1.txt" using 1:2 w l lw 1.5 lc rgb 'blue'
set title "{/:Bold Upwind Method to Inviscid Burgers Equation}"
plot "output2.txt" using 1:2 w l lw 1.5 lc rgb 'red'
unset multiplot