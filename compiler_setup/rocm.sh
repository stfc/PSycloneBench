# Build settings for the LLVM compiler
# ================================================
# Fortran compiler
F90=/opt/rocm-5.0.1/llvm/bin/flang
# C and C++ compiler
CC=/opt/rocm-5.0.1/llvm/bin/clang
CXX=/opt/rocm-5.0.1/llvm/bin/clang++

# C and C++ flags
CFLAGS="-O3 -g"
# Fortran compiler flags
F90FLAGS="-O3 -march=native"
# Flags to use when compiling with OpenMP support
OMPFLAGS=""
# Flags to use when compiling with OpenMP support
#OMPTARGETFLAGS="–fopenmp-targets=nvptx64-nvidia-cuda"
OMPTARGETFLAGS="-O3 -target x86_64-pc-linux-gnu -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx906"

# Linker flags
LDFLAGS="-O3 -target x86_64-pc-linux-gnu -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx906 -lm"
# Location of various CUDA maths libraries
#LDFLAGS+=" -L${CUDA_MATH_DIR}/lib64"
LDFLAGS+=""

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
