# initial config
set term postscript eps enhanced
set output 'roofline.eps'
#set term pngcairo
#set output 'roofline.png'

set nokey
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Operational Intensity (Flops/byte)"
set ylabel "GFlops/s"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

# label offsets
L_MEM_X=0.125
L_MEM_ANG=36

# range of each axis
MAX_X=8
MIN_Y=0.5
MAX_Y=32
set xrange [0.125:MAX_X]
set yrange [MIN_Y:MAX_Y]

# Kernel constants
# First loop nest of shallow has AI = 0.3 FLOP/byte
SHALLOW_LOOP1_AI = 0.3
# u-momentum kernel of nemolite2d has AI = 0.44 FLOP/byte
NEMOLITE_MOM_AI = 0.44

# CPU CONSTANTS
# For single core of Xeon E5-2697 v2 (Archer), as measured with 
# the Intel MKL version of linpack. This is therefore using
# 256-bit AVX instructions (SIMD)
PEAK_GFLOPS=24.1
NUM_CORES=1

#ceilings
C_ALL_CORES		= 1
C_MUL_ADD_BAL	= NUM_CORES
# For Ivy Bridge, AVX registers are 256-bit and therefore can
# hold 4*64-bit double-precision reals. We therefore assume
# that peak, non-SIMD performance is 1/4 that of the performance
# obtained by Linpack
C_SIMD			= 4.0
C_ILP_ONLY		= 2 * C_SIMD

# MEM CONSTANTS
# For single core of Xeon E5-2697 v2 (Archer) as measured with 
# the 'copy' result of STREAM
# with arrays of 15M elements. Therefore, this is bandwidth to 
# main memory, not cache. Units are GB/s.
PEAK_MEM_BW=8.4
# Using arrays of 0.5M elements I think we get bandwidth to
# L3 cache:
PEAK_L3_BW=17.7


NUM_CHANNELS=2
# first ceiling, without multiple memory channels
C_NO_MULTI_CHANNEL	= NUM_CHANNELS

# FUNCTIONS
mem_roof(x,peak)= x * peak
cpu_roof	= PEAK_GFLOPS
min(x, y)	= (x < y) ? x : y
max(x, y)       = (x > y) ? x : y

PEAK_BW = max(PEAK_MEM_BW,PEAK_L3_BW)

cpu_ceiling(x, y)	= min(mem_roof(x,PEAK_BW), y)
mem_ceiling(x)		= min(x, PEAK_GFLOPS)
roofline(x, y)		= cpu_ceiling(x, y)


# LINE STYLES
LINE_ROOF=1
LINE_CEIL=2
LINE_LOOP1_512=3
LINE_LOOP1_1024=4
LINE_MOM_512=5
LINE_MOM_256=6
LINE_MOM_128=7

# Width of the vertical 'bars' at x=1
BAR_WIDTH = 12

set style line LINE_ROOF	lt 1 lw 6 lc rgb "#8B0000"
set style line LINE_CEIL	lt 1 lw 3 lc rgb "blue"

set style line LINE_LOOP1_512     lt 1 lc rgb "dark-olivegreen"
set style line LINE_LOOP1_1024    lt 1 lc rgb "green"
set style line LINE_MOM_512       lt 1 lc rgb "violet"
set style line LINE_MOM_256       lt 1 lc rgb "orange"
set style line LINE_MOM_128       lt 1 lc rgb "red"

# PLOTS
set multiplot

# Bars for measured individual kernel performance

# From Shallow with the Cray compiler (as that's the best)

# Loop1 of shallow with 512^2 achieves 7.0 GFLOPS
set label 12 "shallow: loop 1, 512" at (SHALLOW_LOOP1_AI*0.6),8.0 front textcolor ls LINE_LOOP1_512
set arrow from SHALLOW_LOOP1_AI,MIN_Y to SHALLOW_LOOP1_AI,7.0 nohead ls LINE_LOOP1_512 lw BAR_WIDTH*SHALLOW_LOOP1_AI

set label 13 "shallow: loop 1, 1024" at (SHALLOW_LOOP1_AI*1.06), 4.3 front textcolor ls LINE_LOOP1_1024
# Loop1 of shallow with 1024^2 achieves 4.1 GFLOPS
set arrow from SHALLOW_LOOP1_AI,MIN_Y to SHALLOW_LOOP1_AI,4.1 nohead ls LINE_LOOP1_1024 lw BAR_WIDTH*SHALLOW_LOOP1_AI

# From Nemolite2D with Intel compiler (as that's the fastest)

# 256 domain should fit within L3 cache
set label 14 "nemolite2d: Mom, 256" at (NEMOLITE_MOM_AI*1.06),3.6 front textcolor ls LINE_MOM_256
set arrow from NEMOLITE_MOM_AI,MIN_Y to NEMOLITE_MOM_AI,3.6 nohead ls LINE_MOM_256 lw BAR_WIDTH*NEMOLITE_MOM_AI
# 128 domain - not as fast as you'd expect
set label 15 "nemolite2d: Mom, 128" at (NEMOLITE_MOM_AI*1.06),3.15 front textcolor ls LINE_MOM_128
set arrow from NEMOLITE_MOM_AI,MIN_Y to NEMOLITE_MOM_AI,3.39 nohead ls LINE_MOM_128 lw BAR_WIDTH*NEMOLITE_MOM_AI
# 512 domain ~spills from L3 cache to main memory
set label 11 "nemolite2d: Mom, 512" at (NEMOLITE_MOM_AI*1.06),2.7 front textcolor ls LINE_MOM_512
set arrow from NEMOLITE_MOM_AI,MIN_Y to NEMOLITE_MOM_AI,3.26 nohead ls LINE_MOM_512 lw BAR_WIDTH*NEMOLITE_MOM_AI


# CPU CEILINGS
# All cores (same as roofline)
#set label 3 "All cores used" at (MAX_X-1),(cpu_roof/1.1) right
#plot cpu_ceiling(x, cpu_roof / C_ALL_CORES) ls LINE_CEIL

# SIMD
set label 5 "No SIMD" at (MAX_X-1),((cpu_roof / C_SIMD)/1.1) right
plot cpu_ceiling(x, cpu_roof / C_SIMD) ls LINE_CEIL

# No parallelism
#set label 6 "ILP Only" at (MAX_X-1),((cpu_roof / C_ILP_ONLY)/1.1) right
#plot cpu_ceiling(x, cpu_roof / C_ILP_ONLY) ls LINE_CEIL

# MEM CEILINGS

set label 8 "Main memory" at (L_MEM_X),(mem_roof(L_MEM_X,PEAK_MEM_BW)*1.1) rotate by L_MEM_ANG
plot mem_ceiling(mem_roof(x,PEAK_MEM_BW)) ls LINE_CEIL

# ROOFLINE
set label 1 "Peak FP Performance" at (MAX_X-1),(PEAK_GFLOPS*1.1) right
set label 2 "L3 Mem Bandwidth" at L_MEM_X,mem_roof(L_MEM_X,PEAK_BW)*1.1 rotate by L_MEM_ANG
plot roofline(x, cpu_roof) ls LINE_ROOF

unset multiplot
