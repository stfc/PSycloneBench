# Makefile for the version of the tra_adv benchmark that has had
# the compute moved into multiple, separate subroutines for the
# purposes of processing with Dusk/Dawn.

# Compiler and associated flags must be set through the following
# environment variables, e.g.:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp
#   export MPIF90=mpif90
#
# Scripts are provided to do this in PSycloneBench/compiler_setup

.PHONY: tra_adv_serial tra_adv_acc tra_adv_sir \
        tra_adv_sir_acc

all: tra_adv_serial

# Location of the dl_timer code
DL_TIMER_DIR = ../../../../shared/dl_timer
DL_TIMER_NAME = libdl_timer_omp.a

# Shorthand for invoking PSyclone.
PSYCLONE = psyclone -api nemo -l output

# Serial version.
tra_adv_serial: dl_timer
	mkdir -p $@
	cp tra_adv_driver.F90 $@/.
	cp tra_adv_compute.F90 $@/tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# OpenACC version with timer around outer loop only.
tra_adv_acc: dl_timer
	mkdir -p $@
	${PSYCLONE} -s ../scripts/acc_kernels_unified_memory_trans.py -opsy \
            $@/tra_adv_compute.f90 ./tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	cp tra_adv_driver.F90 $@/.
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} ${UMEMFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ${UMEMFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# OpenACC version with nvtx profiling instrumentation.
tra_adv_acc_prof: dl_timer
ifndef PSYCLONE_NVIDIA_LIB_DIR
	$(error The tra_adv_acc_prof must have the PSYCLONE_NVIDIA_LIB_DIR defined)
endif
	mkdir -p $@
	${PSYCLONE} --profile invokes -s ../scripts/kernels_trans.py -opsy \
            $@/tra_adv_compute.f90 ./tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	cp tra_adv_driver.F90 $@/.
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} ${UMEMFLAGS} -I${PSYCLONE_NVIDIA_LIB_DIR} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ${UMEMFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME} ${PSYCLONE_NVIDIA_LIB_DIR}/libnvtx_prof.a" -C $@

# Serial Fortran version after transformation to SIR-compliant form.
tra_adv_sir: dl_timer
	mkdir -p $@
	${PSYCLONE} -s ../scripts/sir_trans.py -opsy $@/tra_adv_compute.f90 \
            ./tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	cp tra_adv_driver.F90 $@/.
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# OpenACC added after transformation to SIR-compliant form.
tra_adv_sir_acc: dl_timer
	mkdir -p $@
	${PSYCLONE} -s ../scripts/sir_kernels_trans.py -opsy \
            $@/tra_adv_compute.f90 ./tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	cp tra_adv_driver.F90 $@/.
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} ${UMEMFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ${UMEMFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

dl_timer:
	${MAKE} -C ${DL_TIMER_DIR} sm_lib

clean:
	-${MAKE} -C tra_adv_serial $@

allclean: clean
	${MAKE} -C ${DL_TIMER_DIR} allclean
	rm -rf tra_adv_serial
	rm -rf tra_adv_acc
	rm -rf tra_adv_acc_prof
	rm -rf tra_adv_sir
	rm -rf tra_adv_sir_acc
