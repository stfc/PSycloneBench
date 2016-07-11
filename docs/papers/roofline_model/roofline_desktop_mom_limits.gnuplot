# initial config
set term postscript eps enhanced color
set output 'roofline_desktop_mom_limits.eps'

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
# u-momentum kernel of nemolite2d has AI = 0.42 FLOP/byte
NEMOLITE_MOM_AI = 0.42

# CPU CONSTANTS
# For single core of Xeon E5-1620 v2 (my desktop), as measured with 
# the Intel MKL version of linpack. This is therefore using
# 256-bit AVX instructions (SIMD)
PEAK_GFLOPS=28.32

# Upper and lower bounds on performance of u-momentum kernel as
# obtained by analysing the DAG.
C_UMOM_PERFECT_ILP = 2.78
C_UMOM_NO_ILP = 1.73

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
LINE_MOM=5
LINE_MOM_NO_VEC=6
LINE_MOM_CEIL=16
LINE_MOM_SSE_CEIL=17

# Width of the bars
BAR_WIDTH = 0.02

set style line LINE_ROOF	lt 1 lw 6 lc rgb "#8B0000"
set style line LINE_CEIL	lt 1 lw 3 lc rgb "blue"

MOM_COL = "dark-green"
MOM_NO_VEC_COL = "green"
CONT_COL = "dark-blue"
CONT_NO_VEC_COL = "blue"

set style line LINE_MOM       lt 1 lc rgb MOM_COL
set style line LINE_MOM_NO_VEC lt 1 lc rgb MOM_NO_VEC_COL
set style line LINE_MOM_CEIL       lt 1 lw 2 lc rgb MOM_NO_VEC_COL
set style line LINE_MOM_SSE_CEIL   lt 1 lw 2 lc rgb MOM_COL

# PLOTS
set multiplot

# Bars for measured individual kernel performance (GFLOPS)

###########################################################################
# u-Momentum kernel from Nemolite2D with Intel compiler (as that's the fastest)

# 256 domain with SSE
set label 1 "u-Momentum, SSE" at (NEMOLITE_MOM_AI*1.06),2.87 front textcolor ls LINE_MOM
set object 2 rect from (1.0-BAR_WIDTH)*NEMOLITE_MOM_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_MOM_AI,2.87 back fc rgb MOM_COL fs solid

# 256 domain without SSE
set label 2 "u-Momentum" at (NEMOLITE_MOM_AI*1.06),2.09 front textcolor ls LINE_MOM_NO_VEC
set object 3 rect from (1.0-BAR_WIDTH)*NEMOLITE_MOM_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_MOM_AI,2.09 back fc rgb MOM_NO_VEC_COL fs solid

# CPU CEILINGS

# ILP and SIMD

# u-momentum upper bound (perfect ILP)
set label 20 "Perfect ILP, no SIMD" at (MAX_X-1),(C_UMOM_PERFECT_ILP/1.1) right textcolor ls LINE_MOM_CEIL
plot cpu_ceiling(x, C_UMOM_PERFECT_ILP) ls LINE_MOM_CEIL

# u-momentum lower bound (No ILP)
set label 21 "No ILP, no SIMD" at (MAX_X-1),(C_UMOM_NO_ILP/1.1) right textcolor ls LINE_MOM_CEIL
plot cpu_ceiling(x, C_UMOM_NO_ILP) ls LINE_MOM_CEIL

# u-momentum upper bound (perfect ILP) + perfect SSE
set label 24 "Perfect ILP + Perfect SSE" at (MAX_X-1),(2.0*C_UMOM_PERFECT_ILP/1.1) right textcolor ls LINE_MOM_SSE_CEIL
plot cpu_ceiling(x, 2.0*C_UMOM_PERFECT_ILP) ls LINE_MOM_SSE_CEIL

# u-momentum lower bound (No ILP) + perfect SSE
set label 25 "No ILP  + Perfect SSE" at (MAX_X-1),(2.0*C_UMOM_NO_ILP/1.1) right textcolor ls LINE_MOM_SSE_CEIL
plot cpu_ceiling(x, 2.0*C_UMOM_NO_ILP) ls LINE_MOM_SSE_CEIL

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
