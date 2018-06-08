# Build settings for PGI compiler
# ===============================
# Fortran compiler
F90=pgf90
# C compiler
CC=pgcc
# C compiler flags
CFLAGS="-O3"
# Fortran compiler flags
#F90FLAGS+=" -fcheck=all -fbacktrace -ffpe-trap=invalid -g -O0"
F90FLAGS="-O3"
# Linker flags
LDFLAGS="-O3"
# Flags to use when compiling with OpenMP support
OMPFLAGS="-mp"
# Command to use to create archive of object files
AR=ar
# ==============================
export F90
export F90FLAGS
export OMPFLAGS
export AR
