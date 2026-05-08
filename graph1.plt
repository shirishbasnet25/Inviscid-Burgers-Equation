set grid 
set key off 
set size ratio 0.75
set title "{/:Bold Lax Wendroff Method to Inviscid Burgers Equation}"
set xlabel ("x (m)")
set ylabel ("y (m)")
plot "output1.txt" using 1:2 w lp lw 2