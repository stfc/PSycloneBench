# Makefile for the original form of the tra_adv benchmark.

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

tra_adv.exe: timer_lib tra_adv.o
	${F90} ${OMPFLAGS} -o $@ tra_adv.o ${TIMER_DIR}/${TIMER_NAME}

%.o: %.f90
	${F90} ${FORT_FLAGS} -c $<
%.o: %.F90
	${F90} ${FORT_FLAGS} -c $<

timer_lib:
	make -C ${TIMER_DIR} sm_lib

clean:
	rm -f *.o
	rm -f *~
	rm -f tra_adv.exe
