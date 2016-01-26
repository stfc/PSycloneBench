# initial config
set term postscript eps enhanced color
set output 'stream2_desktop.eps'

set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Array size (no. of elements)"
set ylabel "GB/s"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

#set nokey

# range of each axis
#set format x "%e"
set xtics(30,100,500,1000,10000,1000000)
set ytics(142,111,64.0,46.7,20.500,22.5)
set yrange[20.0:150.0]
set multiplot
dfile1='stream2_desktop_fast.dat'
dfile2='stream2_desktop_no_ldfast.dat'
plot dfile1 u 1:10:9:11 w yerrorbars notitle, dfile1 u 1:10 w l lc rgb "red" title "-O3 -fast", dfile2 u 1:10:9:11 w yerrorbars notitle, dfile2 u 1:10 w l lc rgb "blue" title "No -fast"
#dfile='stream2_desktop_no_fast.dat'
#plot dfile u 1:10:9:11 w yerrorbars, dfile u 1:10 w l lc rgb "green"
