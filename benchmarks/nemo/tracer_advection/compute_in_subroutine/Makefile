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

all: tra_adv_serial

# Location of the dl_timer code
DL_TIMER_DIR = ../../../../shared/dl_timer
DL_TIMER_NAME = libdl_timer_omp.a

tra_adv_serial: dl_timer
	mkdir -p $@
	cp tra_adv_*.F90 $@/.
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

dl_timer:
	${MAKE} -C ${DL_TIMER_DIR} sm_lib

clean:
	${MAKE} -C tra_adv_serial $@

allclean: clean
	rm -rf tra_adv_serial
