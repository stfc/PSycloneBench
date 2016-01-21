# initial config
set term postscript eps enhanced color
set output 'stream2_desktop.eps'

set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Array size (no. of elements)"
set ylabel "MB/s"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

set nokey

# range of each axis
#set format x "%e"
set xtics(30,100,500,1000,10000,1000000)
set ytics(150000,117000,65800,46700,20500,20000,2000)
dfile='stream2_desktop.dat'
plot dfile u 1:5 w lp