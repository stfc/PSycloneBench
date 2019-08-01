# Build settings for gfortran compiler
F90=gfortran
CC=gcc

CFLAGS="-O3"
F90FLAGS="-Wall -Wsurprising -Wuninitialized"
#F90FLAGS += -O0
#F90FLAGS += -fcheck=all -fbacktrace -ffpe-trap=invalid -g
F90FLAGS+=" -faggressive-function-elimination"
F90FLAGS+=" -Ofast -mtune=native -finline-limit=50000 -fopt-info-all=gnu_opt_report.txt"
F90FLAGS+=" -march=core2 -mtune=core2"
# f2py does not break long lines so tell gfortran not to
# limit the length of a line
F90FLAGS+=" -ffree-line-length-none"
#F90FLAGS = -O3

OMPFLAGS=""
OMPFLAGS+=" -fopenmp"

LDFLAGS=""

AR=ar

export F90
export F90FLAGS
export OMPFLAGS
export LDFLAGS
export AR
export CC
