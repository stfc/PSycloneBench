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

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a

COMMON_MODULES = clfortran.o opencl_utils_mod.o kernel_args_mod.o \
                 gocean2d_io_mod.o \
                 model_mod.o initialisation_mod.o nemolite2d.o ${INF_LIB}

timer_lib:
	make -C ${TIMER_DIR} sm_lib

inf_lib:
	make -C ${INF_DIR}

nemolite2d.exe: inf_lib ${COMMON_MODULES}
	${F90} -o $@ ${COMMON_MODULES} ${OPENCL_LIBS} ${TIMER_LIB} $(LDFLAGS)

%.o: %.f90
	${F90} ${F90FLAGS} -c -I${INF_INC} -I${TIMER_INC} $<

clean:
	rm -f *.o *~

allclean: clean
	rm -f *.exe *.mod
