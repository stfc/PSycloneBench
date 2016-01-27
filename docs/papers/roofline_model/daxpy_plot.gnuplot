set term postscript eps enhanced color
set output 'daxpy_plot.eps'

set xlabel "Array size (no. of elements)"
set ylabel "GB/s"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

set xrange [30:4000000]
set yrange [20:165]

set grid

set xtics   (30.0000, 100.000, 1000.00, 10000.0, 100000.0, 1.00000e+06)
set ytics   (22.0000, 45.0000, 63.0000, 160.000)

plot 'daxpy_desktop.dat' u 1:4 w lp title 'daxpy', 'daxpypxy_desktop.dat' u 1:4 w lp title 'daxpypxy', 'daxpypxyy_desktop.dat' u 1:4 w lp title 'daxpypxyy', 'daxpypxyyy_destop.dat' u 1:4 w lp title 'daxpypyyy'


