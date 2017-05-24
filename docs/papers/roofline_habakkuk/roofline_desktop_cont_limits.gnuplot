# initial config
set term postscript eps enhanced color
set output 'roofline_desktop_cont_limits.eps'

set nokey
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Operational Intensity (FLOPs/byte)"
set ylabel "GFLOPS"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

# label offsets
L_MEM_X=0.2
L_MEM_ANG=34

# range of each axis
MAX_X=4
MIN_Y=0.5
MAX_Y=34
set xrange [0.1:MAX_X]
set yrange [MIN_Y:MAX_Y]

# Kernel constants
# The continuity kernel of nemolite2d
NEMOLITE_CONT_AI = 0.1276

# CPU CONSTANTS
# For single core of Xeon E5-1620 v2 (my desktop), as measured with 
# the Intel MKL version of linpack. This is therefore using
# 256-bit AVX instructions (SIMD)
PEAK_GFLOPS=28.32

# Upper and lower bounds on performance of u-momentum kernel as
# obtained by analysing the DAG.
C_UMOM_PERFECT_ILP = 2.78
C_UMOM_NO_ILP = 1.73
# Ditto for Continuity kernel
C_CONT_PERFECT_ILP = 4.15
C_CONT_NO_ILP = 2.57

#ceilings
C_ALL_CORES		= 1
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

# FUNCTIONS
mem_roof(x,peak)= x * peak
cpu_roof	= PEAK_GFLOPS
min(x, y)	= (x < y) ? x : y
max(x, y)       = (x > y) ? x : y

PEAK_BW = max(PEAK_MEM_BW,PEAK_L2_BW)

cpu_ceiling(x, y)	= min(mem_roof(x,PEAK_BW), y)
mem_ceiling(x)		= min(x, PEAK_GFLOPS)
roofline(x, y)		= cpu_ceiling(x, y)


# LINE STYLES
LINE_ROOF=1
LINE_CEIL=2
LINE_CONT=10
LINE_CONT_NO_VEC=11
LINE_CONT_CEIL=18
LINE_CONT_SSE_CEIL=19

# Width of the bars
BAR_WIDTH = 0.02

set style line LINE_ROOF	lt 1 lw 6 lc rgb "#8B0000"
set style line LINE_CEIL	lt 1 lw 3 lc rgb "blue"

MOM_COL = "dark-green"
MOM_NO_VEC_COL = "green"
CONT_COL = "dark-blue"
CONT_NO_VEC_COL = "blue"

set style line LINE_CONT        lt 1 lc rgb CONT_COL
set style line LINE_CONT_NO_VEC lt 1 lc rgb CONT_NO_VEC_COL
set style line LINE_CONT_CEIL       lt 1 lw 2 lc rgb CONT_NO_VEC_COL
set style line LINE_CONT_SSE_CEIL   lt 1 lw 2 lc rgb CONT_COL

# PLOTS
set multiplot

# Bars for measured individual kernel performance (GFLOPS)

###########################################################################
# Nemolite2d, Continuity kernel

# 256 domain, SSE
set label 3 "Continuity, SSE" at (NEMOLITE_CONT_AI*1.06),6.05 front textcolor ls LINE_CONT
set object 4 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,6.05 back fc rgb CONT_COL fs solid

# 256 domain, no-vec: 4.13 GFLOPS on desktop
set label 4 "Continuity" at (NEMOLITE_CONT_AI*1.06),4.13 front textcolor ls LINE_CONT_NO_VEC
set object 5 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,4.13 back fc rgb CONT_NO_VEC_COL fs solid

# CPU CEILINGS

# ILP and SIMD

# Continuity upper bound (perfect ILP)
set label 22 "Perfect ILP, no SIMD" at (MAX_X-1),(C_CONT_PERFECT_ILP/1.1) right textcolor ls LINE_CONT_CEIL
plot cpu_ceiling(x, C_CONT_PERFECT_ILP) ls LINE_CONT_CEIL

# Continuity lower bound (No ILP)
set label 23 "No ILP, no SIMD" at (MAX_X-1),(C_CONT_NO_ILP/1.1) right textcolor ls LINE_CONT_CEIL
plot cpu_ceiling(x, C_CONT_NO_ILP) ls LINE_CONT_CEIL

# Continuity lower bound (No ILP) + perfect SSE
set label 26 "No ILP + Perfect SSE" at (MAX_X-1),(2.0*C_CONT_NO_ILP/1.1) right textcolor ls LINE_CONT_SSE_CEIL
plot cpu_ceiling(x, 2.0*C_CONT_NO_ILP) ls LINE_CONT_SSE_CEIL

# Continuity: perfect ILP + perfect SSE
set label 27 "Perfect ILP + Perfect SSE" at (MAX_X-1),(2.0*C_CONT_PERFECT_ILP/1.1) right textcolor ls LINE_CONT_SSE_CEIL
plot cpu_ceiling(x, 2.0*C_CONT_PERFECT_ILP) ls LINE_CONT_SSE_CEIL

# MEM CEILINGS

set label 13 "Main memory" at (L_MEM_X),(mem_roof(L_MEM_X,PEAK_MEM_BW)*1.1) rotate by L_MEM_ANG
plot mem_ceiling(mem_roof(x,PEAK_MEM_BW)) ls LINE_CEIL

set label 17 "L3 cache" at (L_MEM_X),(mem_roof(L_MEM_X,PEAK_L3_BW)*1.1) rotate by L_MEM_ANG
plot mem_ceiling(mem_roof(x,PEAK_L3_BW)) ls LINE_CEIL

# ROOFLINE
set label 14 "Peak FP Performance" at (MAX_X-1),(PEAK_GFLOPS*1.1) right
set label 15 "L2 cache" at L_MEM_X,mem_roof(L_MEM_X,PEAK_BW)*1.1 rotate by L_MEM_ANG
plot roofline(x, cpu_roof) ls LINE_ROOF

unset multiplot
