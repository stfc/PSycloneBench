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
# managed memory
F90FLAGS+=" -acc=gpu -gpu=mem:managed"
# Linker flags
# managed memory
LDFLAGS="-acc=gpu -gpu=mem:managed"
# Location of various CUDA maths libraries. nvtx3interop is required when
# using nvtx for profiling.
LDFLAGS+=" -cuda -L${CUDA_MATH_DIR}/lib64 -lnvtx3interop"
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

