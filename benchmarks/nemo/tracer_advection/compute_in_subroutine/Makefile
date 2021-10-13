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

all: tra_adv.exe

# Location of the dl_timer code
TIMER_DIR = ../../../../shared/dl_timer
TIMER_NAME = libdl_timer_omp.a

FORT_FLAGS = ${F90FLAGS} ${OMPFLAGS} -I${TIMER_DIR}/src

tra_adv.exe: timer_lib tra_adv_compute.o tra_adv_driver.o
	${F90} ${OMPFLAGS} -o $@ tra_adv_compute.o tra_adv_driver.o ${TIMER_DIR}/${TIMER_NAME}

%.o: %.f90
	${F90} ${FORT_FLAGS} -c $<
%.o: %.F90
	${F90} ${FORT_FLAGS} -c $<

timer_lib:
	make -C ${TIMER_DIR} sm_lib

clean:
	rm -f *.o *.mod
	rm -f *~
	rm -f tra_adv.exe

allclean: clean
	rm -f gnu_opt_report.txt output.dat
