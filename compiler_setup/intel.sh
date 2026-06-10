# Build settings for the Intel compiler
F90=ifx
CC=icx

CFLAGS="-O3 -xHost -qopt-report"

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

# For output of source-code-annotated assembler and reports
F90FLAGS+=" -qopt-report=3 -qopt-report-phase=loop,vec"

# Flags to switch-on OpenMP support in compiler
OMPFLAGS="-fopenmp"

LDFLAGS= 
#LDFLAGS+= -fast

# The archiver used to generate the API library.
AR=ar
ARFLAGS=cru

export F90
export F90FLAGS
export CC
export CFLAGS
export OMPFLAGS
export AR
export ARFLAGS
