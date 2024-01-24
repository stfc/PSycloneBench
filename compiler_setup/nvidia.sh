# Build settings for the Nvidia compiler
# ================================================
# Fortran compiler
F90=nvfortran
# C and C++ compiler
CC=nvc
CXX=nvc++

# C and C++ flags
CXXFLAGS="-O3 -Minfo=all"
CFLAGS="-O3"
# Fortran compiler flags
F90FLAGS="-O3 -Minfo=all"

# Debugging options
# The below flags used alongside the above
# would provide debugging information if required
# -Minfo=all provides compiler feedback
# -g provides debugging info to use with GDB
# -O0 tells the compiler to make no optimisations
# -fcheck=all enables run-time tests
# -fbacktrace provides a backtrace of encountered errors
# -ffpe-trap=invalid indicates the line where an error occurs

# Flag to use when compiling with OpenMP support
OMPFLAGS="-mp"
# Flag to use when compiling with OpenMP GPU offloading support
OMPTARGETFLAGS="-mp=gpu -gpu=ccnative"
# Flag to use to specify use of 'managed memory' (unified memory)
UMEMFLAGS="-gpu=managed"
# Flags to use when compiling with OpenACC support
ACCFLAGS="-acc=gpu -gpu=ccnative"

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
export UMEMFLAGS
export ACCFLAGS

export CXXFLAGS
export CFLAGS
export F90FLAGS

export LDFLAGS
export AR
