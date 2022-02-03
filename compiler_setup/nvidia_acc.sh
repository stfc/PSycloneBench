# Build settings for the NVIDIA compiler with OpenACC
# ===================================================
# Fortran compiler
F90=nvfortran
# C compiler
CC=nvcc
# Gives more debug information without performance cost
CFLAGS="-g"
# Fortran compiler flags. -Minfo=all gets the compiler to report on all
# optimisation/parallelisation that it performs.
F90FLAGS="-O3 -Minfo=all"
# Debugging options
#F90FLAGS"+=" -fcheck=all -fbacktrace -ffpe-trap=invalid -g -O0"
# -Mcuda is for CUDA Fortran
# nordc - do not link to routines compiled for device (ensure
# kernel code is in-lined in loops)
# cc = compute capability
# Registers are shared by threads in an SMP. The more registers a kernel
# uses, the fewer threads it can support. This parameter can be tuned and
# should be a multiple of 8.
# -Mcuda is required to build CUDA Fortran
# For Quadro K600
#F90FLAGS+=" -acc -ta=tesla:cc30,nordc -Mcuda=cc30,nordc"
# For Tesla K20c
#F90FLAGS+=" -acc -ta=tesla,cc35,maxregcount:80,nordc -Mcuda=cc35,maxregcount:80,nordc"
# V100 with managed memory
F90FLAGS+=" -acc=gpu -gpu=cc70,managed"
# Linker flags
# For Quadro K600
#LDFLAGS+=" -acc -ta=tesla,cc30 -Mcuda=cc30,nordc"
# For Tesla K20c
#LDFLAGS="-acc -ta=nvidia,cc35 -Mcuda=cc35,nordc"
# V100 with managed memory
LDFLAGS="-acc=gpu -gpu=cc70,managed"
# Location of various CUDA maths libraries. libnvToolsExt is required when
# using nvtx for profiling.
LDFLAGS+=" -Mcuda -L${CUDA_MATH_DIR}/lib64  -lnvToolsExt"
# Flags to use when compiling with OpenMP support
OMPFLAGS="-mp"
# Command to use to create archive of object files
AR=ar
# Location of PSyclone NVIDIA profiling library (used when adding
# profiling to the nemo/tracer_advection benchmark).
PSYCLONE_NVIDIA_LIB_DIR=/home/aporter/PSyclone/lib/profiling/nvidia
# ==============================
export F90
export F90FLAGS
export LDFLAGS
export AR
export PSYCLONE_NVIDIA_LIB_DIR

