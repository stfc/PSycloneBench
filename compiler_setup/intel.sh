# Build settings for the Intel compiler
F90=ifort
CC=icc

CFLAGS="-O3 -axCORE-AVX2 -opt-report"

F90FLAGS=
# -fp-model strict
#F90FLAGS+=" -g -check all -traceback"
#F90FLAGS+=" -O0"
#F90FLAGS+=" -O1"
F90FLAGS+=" -O3"
#F90FLAGS+=" -O4"
#-fast

# SIMD vectorisation and alignment
#F90FLAGS+=" -align array64byte"
#F90FLAGS+=" -no-vec"
#F90FLAGS+=" -axSSE4.2"
#F90FLAGS+=" -xHost"

# Profiling
#F90FLAGS+=" -g -profile-loops=all"
# -guide=4

# Profile-guided optimisation
#F90FLAGS+=" -prof-gen -prof-dir/tmp/profiled"
#F90FLAGS+=" -prof-use -opt-report-phase=pgo"

# Turn-off all compiler limits regarding in-lining of code
F90FLAGS+=" -no-inline-min-size -no-inline-max-per-compile -no-inline-factor"
#F90FLAGS+=" -fno-inline -fno-inline-functions -no-ipo"

# For output of source-code-annotated assembler and reports
#F90FLAGS+=" -S -fsource-asm -fverbose-asm"
F90FLAGS+=" -qopt-report=5 -qopt-report-phase=loop,vec"

# Flags to switch-on OpenMP support in compiler
OMPFLAGS="-openmp"

LDFLAGS= 
#LDFLAGS+= -fast

# The archiver used to generate the API library. We must
# use Intel's xiar if doing IPO as otherwise the library
# doesn't contain the necessary symbols.
AR=xiar
ARFLAGS=cru

export F90
export F90FLAGS
export CC
export CFLAGS
export OMPFLAGS
export AR
export ARFLAGS
