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
SHARED_DIR = ../../../../../../shared

# The dl_timer library
TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
# The dl_esm_lib infrastructure library
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a
# Shared NEMOLite2D routines
COMMON_DIR = ../../../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a

COMMON_MODULES = clfortran.o opencl_utils_mod.o kernel_args_mod.o \
                 nemolite2d.o ${COMMON_LIB} ${INF_LIB}

nemolite2d.exe: inf_lib ${COMMON_MODULES}
	${F90} -o $@ ${COMMON_MODULES} ${OPENCL_LIBS} ${TIMER_LIB} $(LDFLAGS)

timer_lib:
	make -C ${TIMER_DIR} sm_lib

inf_lib:
	make -C ${INF_DIR}

%.o: %.f90
	${F90} ${F90FLAGS} -c -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} $<

clean:
	rm -f *.o *~

allclean: clean
	rm -f *.exe *.mod
