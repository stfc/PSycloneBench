# Makefile for manual OpenMP version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=gfortran
# export F90FLAGS="-O3"
# export OMPFLAGS=" -fopenmp"
# export AR=ar

F90FLAGS += ${OMPFLAGS}
.PHONY: nemolite2d
.DEFAULT_GOAL := nemolite2d

include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${KERNELS} ${INF_LIB}


# Normal targets
nemolite2d: ${COMMON_LIB} timer_lib inf_lib
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)


# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) $(LDFLAGS) ${OMPFLAGS}


# If we need a .f90 file that doesn't exist then it must be a kernel.
# Create a link to the required file...
%.f90:
	ln -sf ../../kernels/fortran/$@ .

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~
	rm -f gnu_opt_report.txt *.optrpt

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

docs:
	doxygen gocean2d.doxy.config

summary:
	@echo -n "${CC}\t\tcpu\t"
	@OMP_NUM_THREADS=4 ./nemolite2d.exe | awk 'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$4} } \
		END {ORS="\n"; print ""}'
