# Build settings for the PGI compiler with OpenACC
F90=pgf90

F90FLAGS=
LDFLAGS=

#F90FLAGS"+=" -fcheck=all -fbacktrace -ffpe-trap=invalid -g -O0"
# -Mcuda is for CUDA Fortran
# nordc - do not link to routines compiled for device (ensure
# kernel code is in-lined in loops)
# cc = compute capability
# Registers are shared by threads in an SMP. The more registers a kernel
# uses, the fewer threads it can support. This parameter can be tuned and
# shoul be a multiple of 8.
#F90FLAGS+=" -O3 -acc -ta=tesla:cc30,nordc -Minfo=all"
#LDFLAGS+=" -O3 -acc -ta=tesla,cc30 -Mcuda=cc30,nordc"
F90FLAGS+=" -O3 -acc -ta=tesla,cc35,maxregcount:80,nordc -Minfo=all"
LDFLAGS+=" -O3 -acc -ta=nvidia,cc35 -Mcuda=cc35,nordc"

AR=ar

export F90
export F90FLAGS
export LDFLAGS
export AR

