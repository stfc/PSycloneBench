# Build settings for gfortran compiler
export F90=gfortran

#NETCDF_INC = ${HOME}/MyInstalls/GCC/gcc_4.9.3/include
#NETCDF_LIB = ${HOME}/MyInstalls/GCC/gcc_4.9.3/lib64

export TIMER_LIB=${HOME}/Projects/dl_timer/dl_timer_lib.a
export TIMER_INC=${HOME}/Projects/dl_timer/src

F90FLAGS="-Wall -Wsurprising -Wuninitialized"
#F90FLAGS += -O0
#F90FLAGS += -fcheck=all -fbacktrace -ffpe-trap=invalid -g
F90FLAGS+=" -faggressive-function-elimination"
F90FLAGS+=" -Ofast -mtune=native -finline-limit=50000 -fopt-info-all=gnu_opt_report.txt"
F90FLAGS+=" -march=core2 -mtune=core2"
# f2py does not break long lines so tell gfortran not to
# limit the length of a line
F90FLAGS+=" -ffree-line-length-none"
#F90FLAGS = -O3

OMPFLAGS=""
OMPFLAGS+=" -fopenmp"

#F90FLAGS+=-I${NETCDF_INC}
#LDFLAGS+=-L${NETCDF_LIB} -lnetcdff -lnetcdf

AR=ar

export F90FLAGS
export OMPFLAGS
export AR
