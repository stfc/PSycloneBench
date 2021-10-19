# Makefile for the version of the tra_adv benchmark that has had
# the compute moved into a subroutine for the purposes of processing
# with Dusk/Dawn.

# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp
#   export MPIF90=mpif90

.PHONY: tra_adv_serial tra_adv_acc

all: tra_adv_serial

# Location of the dl_timer code
DL_TIMER_DIR = ../../../../shared/dl_timer
DL_TIMER_NAME = libdl_timer_omp.a

# Shorthand for invoking PSyclone.
PSYCLONE = psyclone -api nemo -l output

# Serial version.
tra_adv_serial: dl_timer
	mkdir -p $@
	cp tra_adv_*.F90 $@/.
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# OpenACC version with timer around outer loop only.
tra_adv_acc: dl_timer
	mkdir -p $@
	${PSYCLONE} -s ../scripts/kernels_trans.py -opsy $@/tra_adv_compute.f90 ./tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	cp tra_adv_driver.F90 $@/.
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# OpenACC version with nvtx profiling instrumentation.
tra_adv_acc_prof: dl_timer
	mkdir -p $@
	${PSYCLONE} --profile invokes -s ../scripts/kernels_trans.py -opsy $@/tra_adv_compute.f90 ./tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	cp tra_adv_driver.F90 $@/.
	${MAKE} PROF_LIB_INC="-I${PSYCLONE_NVIDIA_LIB_DIR} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME} ${PSYCLONE_NVIDIA_LIB_DIR}/libnvtx_prof.a" -C $@

# Serial Fortran version after transformation to SIR-compliant form.
tra_adv_sir: dl_timer
	mkdir -p $@
	${PSYCLONE} -s ../scripts/sir_trans.py -opsy $@/tra_adv_compute.f90 ./tra_adv_compute.F90
	cp Makefile_gen $@/Makefile
	cp tra_adv_driver.F90 $@/.
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

dl_timer:
	${MAKE} -C ${DL_TIMER_DIR} sm_lib

clean:
	${MAKE} -C tra_adv_serial $@

allclean: clean
	rm -rf tra_adv_serial
