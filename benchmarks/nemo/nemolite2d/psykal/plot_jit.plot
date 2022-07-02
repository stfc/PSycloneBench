#set title 'Bytes transferred by clSetKernelArg'                       # plot title
set ylabel 'bytes'                          # y-axis label
unset key

set terminal pngcairo dashed size 600,600 noenhanced

set style fill solid
set boxwidth 0.7
set yrange [0:*]
set bmargin 4

set output 'bytes.png'
plot "jit_data.txt" using 0:2:xtic(1) with boxes linecolor rgb "#4682B4"

#set title 'Total time kernel execution'                       # plot title
set ylabel 'time (microseconds)'                          # y-axis label
set output 'time.png'
plot "jit_data.txt" using 0:3:xtic(1) with boxes linecolor rgb "#4682B4", \
     '' using 0:3:4:5 with yerrorbars lc rgb 'black' pt 1 lw 2
