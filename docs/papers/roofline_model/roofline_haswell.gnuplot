# initial config
set term postscript eps enhanced color
set output 'roofline_haswell.eps'

set nokey
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Operational Intensity (FLOPs/byte)"
set ylabel "GFLOPS"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

# label offsets
L_MEM_X=0.3
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
# The continuity kernel of nemolite2d
NEMOLITE_CONT_AI = 0.153

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
C_ILP_ONLY		= 2 * C_SIMD

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
LINE_MOM_256_NO_IF=12
LINE_MOM_128=7
LINE_CONT_64=13
LINE_CONT_128=8
LINE_CONT_128_NO_VEC=9
LINE_CONT_256=10
LINE_CONT_256_NO_VEC=11
LINE_CONT_512_NO_VEC=13

# Width of the bars
BAR_WIDTH = 0.02

set style line LINE_ROOF	lt 1 lw 6 lc rgb "#8B0000"
set style line LINE_CEIL	lt 1 lw 3 lc rgb "blue"

MOM_512_COL         = "violet"
MOM_256_COL         = "orange"
MOM_256_NO_IF_COL   = "dark-red"
MOM_128_COL         = "red"
CONT_64_COL         = "purple"
CONT_128_COL        = "green"
CONT_128_NO_VEC_COL = "dark-chartreuse"
CONT_256_COL        = "dark-khaki"
CONT_256_NO_VEC_COL = "dark-goldenrod"
CONT_512_NO_VEC_COL = "dark-green"

set style line LINE_MOM_512       lt 1 lc rgb MOM_512_COL
set style line LINE_MOM_256       lt 1 lc rgb MOM_256_COL
set style line LINE_MOM_256_NO_IF  lt 1 lc rgb MOM_256_NO_IF_COL
set style line LINE_MOM_128       lt 1 lc rgb MOM_128_COL
set style line LINE_CONT_64       lt 1 lc rgb CONT_64_COL
set style line LINE_CONT_128      lt 1 lc rgb CONT_128_COL
set style line LINE_CONT_128_NO_VEC lt 1 lc rgb CONT_128_NO_VEC_COL
set style line LINE_CONT_256        lt 1 lc rgb CONT_256_COL
set style line LINE_CONT_256_NO_VEC lt 1 lc rgb CONT_256_NO_VEC_COL
set style line LINE_CONT_512_NO_VEC lt 1 lc rgb CONT_512_NO_VEC_COL

# PLOTS
set multiplot

# Bars for measured individual kernel performance (GFLOPS)

###########################################################################
# u-Momentum kernel from Nemolite2D with Intel compiler (as that's the fastest)

# 256 domain run without IF's
KERNEL_AI = 0.4492
set label 1 "u-Momentum, no-IFs, 256" at (KERNEL_AI*1.06),5.8 front textcolor ls LINE_MOM_256_NO_IF
set object 1 rect from (1.0-BAR_WIDTH)*KERNEL_AI,MIN_Y to (1.0+BAR_WIDTH)*KERNEL_AI,5.813 back fc rgb MOM_256_NO_IF_COL fs solid

# 128 domain - not as fast as you'd expect
set label 2 "u-Momentum, 128" at (NEMOLITE_MOM_AI*1.06),3.5 front textcolor ls LINE_MOM_128
set object 2 rect from (1.0-BAR_WIDTH)*NEMOLITE_MOM_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_MOM_AI,3.534 back fc rgb MOM_128_COL fs solid

# 256 domain should fit within L3 cache
set label 3 "u-Momentum, 256" at (NEMOLITE_MOM_AI*1.06),3.6 front textcolor ls LINE_MOM_256
set object 3 rect from (1.0-BAR_WIDTH)*NEMOLITE_MOM_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_MOM_AI,3.737 back fc rgb MOM_256_COL fs solid

# 512 domain ~spills from L3 cache to main memory
set label 4 "u-Momentum, 512" at (NEMOLITE_MOM_AI*1.06),3.0 front textcolor ls LINE_MOM_512
set object 4 rect from (1.0-BAR_WIDTH)*NEMOLITE_MOM_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_MOM_AI,3.504 back fc rgb MOM_512_COL fs solid

###########################################################################
# Nemolite2d, Continuity kernel

# 256 domain, SSE
set label 5 "Continuity, SSE, 256" at (NEMOLITE_CONT_AI*1.06),5.8 front textcolor ls LINE_CONT_256
set object 5 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,5.946 back fc rgb CONT_256_COL fs solid

# 128 domain, SSE
set label 6 "Continuity, SSE, 128" at (NEMOLITE_CONT_AI*1.06),5.6 front textcolor ls LINE_CONT_128
set object 6 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,5.717 back fc rgb CONT_128_COL fs solid

# 64 domain, SSE
set label 7 "Continuity, SSE, 64" at (NEMOLITE_CONT_AI*1.06),5.1 front textcolor ls LINE_CONT_64
set object 7 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,5.251 back fc rgb CONT_64_COL fs solid

# 128 domain, no-vec
#set label 8 "Continuity, no-vec, 128" at (NEMOLITE_CONT_AI*1.06),3.6 front textcolor ls LINE_CONT_128_NO_VEC
#set object 8 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,3.558 back fc rgb CONT_128_NO_VEC_COL fs solid

# 256 domain, no-vec
#set label 9 "Continuity, no-vec, 256" at (NEMOLITE_CONT_AI*1.06),3.2 front textcolor ls LINE_CONT_256_NO_VEC
#set object 9 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,3.410 back fc rgb CONT_256_NO_VEC_COL fs solid

# 512 domain, no-vec
set label 10 "Continuity, no-vec, 512" at (NEMOLITE_CONT_AI*1.06),2.8 front textcolor ls LINE_CONT_512_NO_VEC
set object 10 rect from (1.0-BAR_WIDTH)*NEMOLITE_CONT_AI,MIN_Y to (1.0+BAR_WIDTH)*NEMOLITE_CONT_AI,2.981 back fc rgb CONT_512_NO_VEC_COL fs solid

# CPU CEILINGS
# All cores (same as roofline)
#set label 3 "All cores used" at (MAX_X-1),(cpu_roof/1.1) right
#plot cpu_ceiling(x, cpu_roof / C_ALL_CORES) ls LINE_CEIL

# SIMD
set label 11 "No SIMD" at (MAX_X-1),((cpu_roof / C_SIMD)/1.1) right
plot cpu_ceiling(x, cpu_roof / C_SIMD) ls LINE_CEIL

# No parallelism
#set label 12 "ILP Only" at (MAX_X-1),((cpu_roof / C_ILP_ONLY)/1.1) right
#plot cpu_ceiling(x, cpu_roof / C_ILP_ONLY) ls LINE_CEIL

# MEM CEILINGS

set label 13 "Main memory" at (L_MEM_X),(mem_roof(L_MEM_X,PEAK_MEM_BW)*1.1) rotate by L_MEM_ANG
plot mem_ceiling(mem_roof(x,PEAK_MEM_BW)) ls LINE_CEIL

# ROOFLINE
set label 14 "Peak FP Performance" at (MAX_X-1),(PEAK_GFLOPS*1.1) right
set label 15 "L3 Mem Bandwidth" at L_MEM_X,mem_roof(L_MEM_X,PEAK_BW)*1.1 rotate by L_MEM_ANG
plot roofline(x, cpu_roof) ls LINE_ROOF

unset multiplot
