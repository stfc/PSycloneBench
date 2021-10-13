# Makefile for the original form of the tra_adv benchmark.

# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp
#   export MPIF90=mpif90

.PHONY: all tra_adv_serial

all: tra_adv_serial

# Shorthand for invoking PSyclone
PSYCLONE = psyclone -api nemo -l output

tra_adv_serial:
	mkdir -p $@
	cp tra_adv.F90 $@/.
	cp Makefile_gen $@/Makefile
	${MAKE} -C $@

tra_adv_acc:
	mkdir -p $@
	${PSYCLONE} -s ./kernels_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} -C $@

clean:
	${MAKE} -C tra_adv_serial $@
	${MAKE} -C tra_adv_acc $@

allclean:
	rm -rf tra_adv_serial
	rm -rf tra_adv_acc
