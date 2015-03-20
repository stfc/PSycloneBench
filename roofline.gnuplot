# initial config
set term postscript eps enhanced
set output 'roofline.eps'
set nokey
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set xlabel "Operational Intensity"
set ylabel "GFlops"

# sets log base 2 scale for both axes
set logscale x 2
set logscale y 2

AI_FILE = "shallow_arithmetic_intensity.dat"
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
# Each momentum kernel (u,v) of nemolite2d has AI = 0.44 FLOP/byte
NEMOLITE_MOM_AI = 0.44

# CPU CONSTANTS
# For single core of Xeon E5-2697 v2 (Archer), as measured with 
# the Intel MKL version of linpack
PEAK_GFLOPS=24.1
NUM_CORES=1

#ceilings
C_ALL_CORES		= 1
C_MUL_ADD_BAL	= NUM_CORES
C_SIMD			= 2 * C_MUL_ADD_BAL
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

cpu_ceiling(x, y)	= min(mem_roof(x,PEAK_MEM_BW), y)
mem_ceiling(x)		= min(x, PEAK_GFLOPS)
roofline(x, y)		= cpu_ceiling(x, y)

# LINE STYLES
LINE_ROOF=1
LINE_CEIL=2
LINE_LOOP1_512=3
LINE_LOOP1_1024=4
LINE_MOM_512=5

set style line LINE_ROOF	lt 1 lw 6 lc rgb "#8B0000"
set style line LINE_CEIL	lt 1 lw 3 lc rgb "blue"
set style line LINE_LOOP1_512       lt 1 lw 4 lc rgb "red"
set style line LINE_LOOP1_1024       lt 1 lw 4 lc rgb "green"
set style line LINE_MOM_512       lt 1 lw 4 lc rgb "violet"

# PLOTS
set multiplot


# CPU CEILINGS
# All cores (same as roofline)
set label 3 "All cores used" at (MAX_X-1),(cpu_roof/1.1) right
plot cpu_ceiling(x, cpu_roof / C_ALL_CORES) ls LINE_CEIL

# MUL/ADD balance / only 1 core
set label 4 "Mul/Add balance, No TLP" at (MAX_X-1),((cpu_roof / C_MUL_ADD_BAL)/1.1) right
plot cpu_ceiling(x, cpu_roof / C_MUL_ADD_BAL) ls LINE_CEIL

# SIMD
set label 5 "With SIMD" at (MAX_X-1),((cpu_roof / C_SIMD)/1.1) right
plot cpu_ceiling(x, cpu_roof / C_SIMD) ls LINE_CEIL

# No paralellism
set label 6 "ILP Only" at (MAX_X-1),((cpu_roof / C_ILP_ONLY)/1.1) right
plot cpu_ceiling(x, cpu_roof / C_ILP_ONLY) ls LINE_CEIL

# MEM CEILINGS
# No dual channel
#set label 7 "No Dual Channel" at (L_MEM_X),(mem_roof(L_MEM_X,PEAK_MEM_BW)/C_NO_MULTI_CHANNEL*1.1) rotate by L_MEM_ANG
#plot mem_ceiling(mem_roof(x,PEAK_MEM_BW) / C_NO_MULTI_CHANNEL) ls LINE_CEIL

set label 8 "In L3" at (L_MEM_X),(mem_roof(L_MEM_X,PEAK_L3_BW)*1.1) rotate by L_MEM_ANG
plot mem_ceiling(mem_roof(x,PEAK_L3_BW)) ls LINE_CEIL


#set arrow from SHALLOW_LOOP1_AI,MIN_Y to SHALLOW_LOOP1_AI,mem_roof(0.3,PEAK_MEM_BW) nohead ls LINE_LOOP1
set arrow from SHALLOW_LOOP1_AI,MIN_Y to SHALLOW_LOOP1_AI,7.0 nohead ls LINE_LOOP1_512
set arrow from SHALLOW_LOOP1_AI,MIN_Y to SHALLOW_LOOP1_AI,4.1 nohead ls LINE_LOOP1_1024

# Momentum kernel of nemolite2d has AI = 0.44 FLOP/byte
set label 11 "nemolite2d: Momentum" at (NEMOLITE_MOM_AI),1.0
set arrow from NEMOLITE_MOM_AI,MIN_Y to NEMOLITE_MOM_AI,mem_roof(NEMOLITE_MOM_AI,PEAK_MEM_BW) nohead ls LINE_MOM_512

# ROOFLINE
set label 1 "Peak FP Performance" at (MAX_X-1),(PEAK_GFLOPS*1.1) right
set label 2 "Peak Mem Bandwidth" at L_MEM_X,mem_roof(L_MEM_X,PEAK_MEM_BW)*1.1 rotate by L_MEM_ANG
plot roofline(x, cpu_roof) ls LINE_ROOF

unset multiplot
