# Build settings for the AMDFLANG compiler
# ================================================
# This is an experimental file so other flags may be
# needed for accelerated compilation
# Alternative flags have been provided in the comments
# where they have been found to be useful
# These flags are for ROCM AFAR 22.2.0.

# Fortran compiler
F90=amdflang
# C and C++ compiler
CC=amdclang
CXX=amdclang++

# C and C++ flags
# note that -g is used for debugging information
# as this is an experimental implementation
CFLAGS="-O3 -g"
# Fortran compiler flags
# As above, -g provides debugging information
F90FLAGS="-O3 -g -fsave-optimization-record -ffast-math"
# Flags to use when compiling with OpenMP support
OMPFLAGS="-fopenmp=libomp"
# Flags to use when compiling with OpenMP GPU offloading support
OMPTARGETFLAGS="-fopenmp=libomp --offload-arch=gfx942 -mp -fopenmp-offload-mandatory -fopenmp-force-usm"

UMEMFLAGS=" "

# Linker flags
LDFLAGS="-lomp -lomptarget -fopenmp=libomp -fno-lto --offload-arch=gfx942 -fopenmp-offload-mandatory -fopenmp-force-usm -lflang_rt.hostdevice"
## Location of various CUDA maths libraries
#LDFLAGS+=" -L${CUDA_MATH_DIR}/lib64"

# Command to use to create archive of object files
AR=ar

# ==============================
export F90
export CC
export CXX

export OMPFLAGS
export OMPTARGETFLAGS
export UMEMFLAGS

export CFLAGS
export F90FLAGS

export LDFLAGS
export AR
