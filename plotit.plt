set title "Fit Results"
set xlabel "x"
set ylabel "y"
plot  "xy_pts.dat" with points title "Node Points"
replot "fit.dat" using 1:2 with lines title "Fit"
set output "plot.ps"
set terminal postscript enhanced color landscape
replot
set output "plot.png"
set terminal png
replot

