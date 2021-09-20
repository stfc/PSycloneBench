# Build settings for the PGI compiler with OpenACC
# ================================================
# Fortran compiler
F90=flang
# C and C++ compiler
CC=clang
CXX=clang++

# C and C++ flags
CFLAGS="-O3 -march=native -g"
# Fortran compiler flags
F90FLAGS="-O3 -march=native -g"
# Flags to use when compiling with OpenMP support
OMPFLAGS="-fopenmp"
# Flags to use when compiling with OpenMP support
#OMPTARGETFLAGS="–fopenmp-targets=nvptx64-nvidia-cuda"
OMPTARGETFLAGS="-fopenmp -fopenmp-targets=nvptx64"

# Linker flags
LDFLAGS="-lomp -lomptarget"
# Location of various CUDA maths libraries
#LDFLAGS+=" -L${CUDA_MATH_DIR}/lib64"
LDFLAGS+=" -L/apps/packages/compilers/llvm/13.0.0rc3/lib"

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

