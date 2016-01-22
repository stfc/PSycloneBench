# initial config
set term postscript eps enhanced color
set output 'stream2_archer.eps'

# Data file
dfile='stream2_archer.dat'

set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Array size (no. of elements)"
set ylabel "MB/s"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

set nokey

# range of each axis
set xtics(30,100,500,1000,10000,1000000)
set ytics(130000,108000,54000,37700,17600)

plot dfile u 1:5 w lp