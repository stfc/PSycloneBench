# Build settings for the Nvidia compiler
# ================================================
# Fortran compiler
F90=nvfortran
# C and C++ compiler
CC=nvc
CXX=nvc++

# C and C++ flags
CFLAGS="-O3 -Minfo=all -g"
# Fortran compiler flags
F90FLAGS="-O3 -Minfo=all"

# Debugging options
# F90FLAGS"+=" -fcheck=all -fbacktrace -ffpe-trap=invalid -g -O0"

# Flag to use when compiling with OpenMP support
OMPFLAGS="-mp"
# Flag to use when compiling with OpenMP GPU offloading support
OMPTARGETFLAGS="-mp=gpu"
# Flags to use when compiling with OpenACC support
ACCFLAGS="-acc -ta=tesla"

# Linker flags
LDFLAGS=""
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
export ACCFLAGS

export CFLAGS
export F90FLAGS

export LDFLAGS
export AR
