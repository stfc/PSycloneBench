# Build settings for PGI compiler
F90=pgf90
CC=pgcc

CFLAGS="-O3"

#F90FLAGS+=" -fcheck=all -fbacktrace -ffpe-trap=invalid -g -O0"

F90FLAGS="-O3"
LDFLAGS="-O3"

OMPFLAGS="-fopenmp"

AR=ar

export F90
export F90FLAGS
export OMPFLAGS
export AR
