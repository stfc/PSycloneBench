# initial config
#set term pngcairo
#set output 'stream2_desktop.png'

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
set xtics(20,100,500,1000,10000,1000000)
set ytics(142,111,65.0,46.7,30,20.70)
set yrange[20.0:150.0]
set multiplot
dfile1='stream2_runs_no_HT.out'
dfile2='stream2_runs_no_HT_no_turb_no_step.out'
plot dfile1 u 1:10 w lp notitle, dfile2 u 1:10 w lp title 'No turbo or stepping'

