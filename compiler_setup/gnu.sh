# Build settings for gfortran compiler
F90=gfortran
CC=gcc
CXX=g++

# Common optimization flags for CFLAGS and F90FLAGS
OPTFLAGS=" -Ofast -mtune=native -finline-limit=50000 -fopt-info-all=gnu_opt_report.txt"

CFLAGS=$OPTFLAGS
F90FLAGS="-Wall -Wsurprising -Wuninitialized"
#F90FLAGS += -O0
#F90FLAGS += -fcheck=all -fbacktrace -ffpe-trap=invalid -g
F90FLAGS+=" -faggressive-function-elimination"
F90FLAGS+=$OPTFLAGS
# f2py does not break long lines so tell gfortran not to
# limit the length of a line
F90FLAGS+=" -ffree-line-length-none"

OMPFLAGS=""
OMPFLAGS+=" -fopenmp"

LDFLAGS=""

AR=ar

export F90
export F90FLAGS
export CFLAGS
export OMPFLAGS
export LDFLAGS
export AR
export CC
export CXX
