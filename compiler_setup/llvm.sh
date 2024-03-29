# Build settings for the LLVM compiler
# ================================================
# This is an experimental file so other flags may be
# needed for accelerated compilation
# Alternative flags have been provided in the comments
# where they have been found to be useful

# Fortran compiler
F90=flang
# C and C++ compiler
CC=clang
CXX=clang++

# C and C++ flags
# note that -g is used for debugging information
# as this is an experimental implementation
CFLAGS="-O3 -march=native -g"
# Fortran compiler flags
# As above, -g provides debugging information
F90FLAGS="-O3 -march=native -g"
# Flags to use when compiling with OpenMP support
OMPFLAGS="-fopenmp"
# Flags to use when compiling with OpenMP GPU offloading support
OMPTARGETFLAGS="-fopenmp -fopenmp-targets=nvptx64"
# OMPTARGETFLAGS="–fopenmp-targets=nvptx64-nvidia-cuda" 

# Linker flags
LDFLAGS="-lomp -lomptarget"
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
