# Build settings for the PGI compiler with OpenACC
# ================================================
# Fortran compiler
F90=pgf90
# C compiler
CC=pgcc
CFLAGS="-g"
# Fortran compiler flags
F90FLAGS="-O3 -Minfo=all"
#F90FLAGS"+=" -fcheck=all -fbacktrace -ffpe-trap=invalid -g -O0"
# -Mcuda is for CUDA Fortran
# nordc - do not link to routines compiled for device (ensure
# kernel code is in-lined in loops)
# cc = compute capability
# Registers are shared by threads in an SMP. The more registers a kernel
# uses, the fewer threads it can support. This parameter can be tuned and
# shoul be a multiple of 8.
# -Mcuda is required to build CUDA Fortran
# For Quadro K600
F90FLAGS+=" -acc -ta=tesla:cc30,nordc -Mcuda=cc30,nordc"
# For Tesla K20c
#F90FLAGS+=" -acc -ta=tesla,cc35,maxregcount:80,nordc -Mcuda=cc35,maxregcount:80,nordc"
# Linker flags
# For Quadro K600
LDFLAGS+=" -acc -ta=tesla,cc30 -Mcuda=cc30,nordc"
# For Tesla K20c
#LDFLAGS="-acc -ta=nvidia,cc35 -Mcuda=cc35,nordc"
# Flags to use when compiling with OpenMP support
OMPFLAGS="-mp"
# Command to use to create archive of object files
AR=ar
# ==============================
export F90
export F90FLAGS
export LDFLAGS
export AR

