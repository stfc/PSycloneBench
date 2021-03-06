# Build settings for gfortran compiler
F90=mpifort
CC=mpicc
CXX=mpic++

# Common optimization flags for CFLAGS and F90FLAGS
OPTFLAGS=" -O3 -mtune=native -finline-limit=50000 -fopt-info-all=gnu_opt_report.txt"
OPTFLAGS+=" -g"

CFLAGS=$OPTFLAGS
F90FLAGS="-Wall -Wsurprising -Wuninitialized"
#F90FLAGS += -O0
#F90FLAGS += -fcheck=all -fbacktrace -ffpe-trap=invalid -g
F90FLAGS+=" -faggressive-function-elimination"
F90FLAGS+=$OPTFLAGS

MCXXFLAGS="--Wn,-O3,-mtune=native,-finline-limit=50000,-fopt-info-all=gnu_opt_report.txt"

OMPFLAGS=""
OMPFLAGS=" -fopenmp"

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
