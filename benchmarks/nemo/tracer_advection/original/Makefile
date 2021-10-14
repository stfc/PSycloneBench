# Makefile for the original form of the tra_adv benchmark. Contains
# targets for building the original (serial) version, a PSyclone-processed
# OpenACC version and a PSyclone-processed OpenMP version.

# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp

.PHONY: all tra_adv_serial tra_adv_acc tra_adv_omp

# The default is to build the original, serial CPU version
all: tra_adv_serial

# Location of the dl_timer code
DL_TIMER_DIR = ../../../../shared/dl_timer
DL_TIMER_NAME = libdl_timer_omp.a

# Shorthand for invoking PSyclone.
PSYCLONE = psyclone -api nemo -l output

tra_adv_serial: dl_timer
	mkdir -p $@
	cp tra_adv.F90 $@/.
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# The OpenACC-accelerated version has the existing dl_timer calls but also
# has profiling instrumentation added that uses NVTX.  We profile 'invokes'
# rather than 'kernels' since the whole outer loop gets put into a single
# ACC KERNELS region.
tra_adv_acc:
	mkdir -p $@
	${PSYCLONE} --profile invokes -s ./kernels_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I${PSYCLONE_NVIDIA_LIB_DIR} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME} ${PSYCLONE_NVIDIA_LIB_DIR}/libnvtx_prof.a" -C $@

tra_adv_omp:
	mkdir -p $@
	${PSYCLONE} -s ./omp_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

dl_timer:
	make -C ${DL_TIMER_DIR} sm_lib

clean:
	${MAKE} -C tra_adv_serial $@
	${MAKE} -C tra_adv_acc $@
	${MAKE} -C tra_adv_omp $@

allclean:
	rm -rf tra_adv_serial
	rm -rf tra_adv_acc
