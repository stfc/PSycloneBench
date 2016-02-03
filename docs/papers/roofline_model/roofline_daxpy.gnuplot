# initial config
set term postscript eps enhanced color
set output 'roofline_daxpy.eps'

set nokey
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Operational Intensity (FLOPs/byte)"
set ylabel "GFLOPS"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

# label offsets
L_MEM_X=0.3
L_MEM_ANG=27

# range of each axis
MAX_X=2
MIN_Y=0.5
MAX_Y=34
set xrange [0.1:MAX_X]
set yrange [MIN_Y:MAX_Y]

# CPU CONSTANTS
# For single core of Xeon E5-1620 v2 (my desktop), as measured with 
# the Intel MKL version of linpack. This is therefore using
# 256-bit AVX instructions (SIMD)
PEAK_GFLOPS=28.32
NUM_CORES=1

#ceilings
C_ALL_CORES		= 1
C_MUL_ADD_BAL	= NUM_CORES
# For Ivy Bridge, AVX registers are 256-bit and therefore can
# hold 4*64-bit double-precision reals. We therefore assume
# that peak, non-SIMD performance is 1/4 that of the performance
# obtained by Linpack
C_SIMD			= 4.0

# MEM CONSTANTS
# For single core of Xeon E5-1620 v2 (desktop) as measured with 
# the 'DAXPY' result of STREAM2. Units are GB/s.
PEAK_MEM_BW=20.5
PEAK_L3_BW=46.0
PEAK_L2_BW=61.0
PEAK_L1_BW=160.0


NUM_CHANNELS=2
# first ceiling, without multiple memory channels
C_NO_MULTI_CHANNEL	= NUM_CHANNELS

# FUNCTIONS
mem_roof(x,peak)= x * peak
cpu_roof	= PEAK_GFLOPS
min(x, y)	= (x < y) ? x : y
max(x, y)       = (x > y) ? x : y

PEAK_BW = max(PEAK_MEM_BW,PEAK_L1_BW)

cpu_ceiling(x, y)	= min(mem_roof(x,PEAK_BW), y)
mem_ceiling(x)		= min(x, PEAK_GFLOPS)
roofline(x, y)		= cpu_ceiling(x, y)

LINE_ROOF=1
LINE_CEIL=2
LINE_CPU_CEIL = 3

# Width of the bars
BAR_WIDTH = 0.02

set style line LINE_ROOF	lt 1 lw 6 lc rgb "black"
set style line LINE_CEIL	lt 1 lw 3 lc rgb "blue"
set style line LINE_CPU_CEIL	lt 1 lw 3 lc rgb "dark-blue"

kernels =          "AXPY AXPYPXY AXPYPXYY AXPYPXYYY AXPYPXYYY"
kernel_ai =       "0.125 0.167    0.208     0.25     0.25"
kernel_flops_L3 = "3.65   7.18 8.59  10.26 9.57"
kernel_flops_L2 = "5.08  10.49 12.49 14.4  12.5 "
kernel_flops_L1 = "13.30 21.74 22.70 21.76 15.59"
kernel_xshift = "0.0 0.0 0.0 -0.01 0.01"
colors = "violet orange dark-red red brown pink"
L3_colour = "black"
L2_colour = "red"
L1_colour = "blue"

set multiplot

# Set up the line types
set for [i=1:words(colors)] linetype i lc rgb word(colors, i)

# Draw a rectangle for each data point
# Problem size fits in L1
obj_idx = 0
set for [i=1:words(kernels)-1] object i+obj_idx rect from (1.0-BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),MIN_Y to (1.0+BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),word(kernel_flops_L1, i) back fc rgb L1_colour fs solid
set for [i=words(kernels):words(kernels)] object i+obj_idx rect from (1.0-BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),MIN_Y to (1.0+BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),word(kernel_flops_L1, i) back fc rgb L1_colour fs pattern 1

# Problem size fits in L2
obj_idx = words(kernels)
set for [i=1:words(kernels)-1] object i+obj_idx rect from (1.0-BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),MIN_Y to (1.0+BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),word(kernel_flops_L2, i) back fc rgb L2_colour fs solid
set for [i=words(kernels):words(kernels)] object i+obj_idx rect from (1.0-BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),MIN_Y to (1.0+BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),word(kernel_flops_L2, i) back fc rgb L2_colour fs pattern 1

# Problem size fits in L3
obj_idx = obj_idx + words(kernels)
set for [i=1:words(kernels)-1] object i+obj_idx rect from (1.0-BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),MIN_Y to (1.0+BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),word(kernel_flops_L3, i) back fc rgb L3_colour fs solid
set for [i=words(kernels):words(kernels)] object i+obj_idx rect from (1.0-BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),MIN_Y to (1.0+BAR_WIDTH+word(kernel_xshift,i))*word(kernel_ai, i),word(kernel_flops_L3, i) back fc rgb L3_colour fs pattern 1

# Label each cluster of bars
xshift = 0.02
# Put a white box behind each label
set for [i=1:words(kernels)] object i+20 rect from (1.0-BAR_WIDTH-xshift)*word(kernel_ai,i),MIN_Y*1.3 to (1.0+BAR_WIDTH+xshift)*word(kernel_ai,i),MIN_Y*3.1 back fc rgb "white" fs solid noborder
# The labels themselves
set for [i=1:words(kernels)] label i+20 word(kernels,i) at word(kernel_ai,i),MIN_Y*2.0 centre rotate by 90

# CPU CEILINGS

# SIMD
set label 11 "No SIMD" at (MAX_X-0.5),((cpu_roof / C_SIMD)/1.1) right tc rgb "dark-blue"
plot cpu_ceiling(x, cpu_roof / C_SIMD) ls LINE_CPU_CEIL

# MEM CEILINGS

set label 13 "Memory Bandwidth" at 0.45,(mem_roof(0.45,PEAK_MEM_BW)*0.87) rotate by L_MEM_ANG tc rgb "blue"
set label 16 "L2 Bandwidth" at (L_MEM_X),(mem_roof(L_MEM_X,PEAK_L2_BW)*0.87) rotate by L_MEM_ANG tc rgb "blue"
set label 17 "L3 Bandwidth" at 0.34,(mem_roof(0.34,PEAK_L3_BW)*0.87) rotate by L_MEM_ANG tc rgb "blue"
plot mem_ceiling(mem_roof(x,PEAK_MEM_BW)) ls LINE_CEIL
plot mem_ceiling(mem_roof(x,PEAK_L3_BW)) ls LINE_CEIL
plot mem_ceiling(mem_roof(x,PEAK_L2_BW)) ls LINE_CEIL
# ROOFLINE
set label 14 "Peak FP Performance" at (MAX_X-0.5),(PEAK_GFLOPS*1.1) right
set label 15 "L1 Bandwidth" at 0.125,mem_roof(0.125,PEAK_BW)*1.1 rotate by L_MEM_ANG
plot roofline(x, cpu_roof) ls LINE_ROOF

unset multiplot
