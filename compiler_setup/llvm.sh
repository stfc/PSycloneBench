# Build settings for the LLVM compiler
# ================================================
# This is an experimental file so other flags may be
# needed for accelerated compilation

# Fortran compiler
F90=flang
# If flang is not available or causes compiler errors uncomment gfortran:
# F90=gfortran
# C and C++ compiler
CC=clang
CXX=clang++

# C and C++ flags
CFLAGS="-O3"
# Fortran compiler flags
F90FLAGS="-O3"
# Flags to use when compiling with OpenMP support
OMPFLAGS="-fopenmp"
# Flags to use when compiling with OpenMP GPU offloading support
# For AMD Rocm (march is MI50: fgx906, MI100: gfx908):
# OMPTARGETFLAGS="-target x86_64-pc-linux-gnu -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx908"
# For NVIDIA:
OMPTARGETFLAGS="–fopenmp-targets=nvptx64-nvidia-cuda"

# Linker flags
LDFLAGS="-fopenmp"

# Location of various CUDA maths libraries
LDFLAGS+=" -L${CUDA_MATH_DIR}/lib64"

# Command to use to create archive of object files
AR=ar

# ==============================
export F90
export CC
export CXX

export OMPFLAGS
export OMPTARGETFLAGS

export CFLAGS
export F90FLAGS

export LDFLAGS
export AR
