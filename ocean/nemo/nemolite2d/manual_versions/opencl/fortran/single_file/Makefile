# This Makefile picks up the compiler to use, the location of the
# OpenCL libraries and any flags from environment variables. e.g. to
# use AMD's dispatch library:
# 
# export F90=gfortran
# export CC=gcc
# export OPENCL_LIBS="-L/opt/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL"
# export OPENCL_INCLUDE=-I/opt/AMDAPPSDK-3.0/include
# export CFLAGS="-g -O0 -Wall"
# export LDFLAGS=-lm

# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../../../../../shared

# The dl_timer library
TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
# The dl_esm_lib infrastructure library
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a
# Shared NEMOLite2D routines
COMMON_DIR = ../../../../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a
# Location of the OpenCL kernels
KERNEL_DIR = ../../../../kernels/opencl

COMMON_MODULES = opencl_utils_mod.o kernel_args_mod.o \
                 nemolite2d.o

nemolite2d.exe: kernels ${COMMON_LIB} inf_lib timer_lib ${COMMON_MODULES}
	${F90} ${OMPFLAGS} -o $@ ${COMMON_MODULES} ${COMMON_LIB} \
${INF_LIB} ${TIMER_LIB} ${OPENCL_LIBS} $(LDFLAGS)

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

kernels:
	${MAKE} -C ${KERNEL_DIR} CC=g++

# The clfortran module need only be compiled for its interface definitions.
# It *must not* be linked in to the final executable as it collides with
# the real OpenCL run-time and nothing works!
opencl_utils_mod.o: clfortran.o

${COMMON_LIB}:
	${MAKE} -C ${COMMON_DIR}

%.o: %.f90
	${F90} ${F90FLAGS} -c -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} $<

clean:
	rm -f *.o *~

allclean: clean
	rm -f *.exe *.mod
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean
	${MAKE} -C ${COMMON_DIR} allclean
