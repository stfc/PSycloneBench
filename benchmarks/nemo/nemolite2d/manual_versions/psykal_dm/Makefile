# Makefile for manual MPI version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler and flags from
# environment variables.

# Make sure that F90 has a MPI compiler selected.
ifeq (,$(findstring mpi,$(F90)))
	F90 := mpif90 -fc=${F90}
endif

.PHONY: all nemolite2d nemolite2d-serial
.DEFAULT_GOAL := nemolite2d
include ../common.mk

all: nemolite2d

# The modules that are common to both targets (INF_LIB is an archive that
# must come at the end of the list of objects passed to the linker)
COMMON_MODULES = ${COMMON_LIB} ${KERNELS} time_step_mod.o ${INF_LIB}

# Interdependencies between modules, alphabetical order
nemolite2d.o: $(COMMON_MODULES)

time_step_mod.o: ${COMMON_LIB} inf_lib_parallel


# Normal targets
nemolite2d: timer_lib_parallel inf_lib_parallel
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" MPI=yes nemolite2d.exe

nemolite2d-serial: timer_lib_serial inf_lib_serial
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" MPI=no nemolite2d.exe

# Generic rules
nemolite2d.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) $(TIMER_LIB) $(LDFLAGS) ${OMPFLAGS}

clean: 
	rm -f *.o *.mod *.MOD *~

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	make -C ${TIMER_DIR} allclean
	make -C ${INF_DIR} distclean
	make -C ${COMMON_DIR} allclean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/
	rm -rf *.dat *.txt

docs:
	doxygen gocean2d.doxy.config

summary:
	@echo -n "${CC}\t\tcpu\t"
	@mpirun -n $(shell nproc) ./nemolite2d.exe | awk \
		'BEGIN {ORS="\t"; found_ua=0; found_va=0} \
		{if ($$1 == "ua" && found_ua == 0) {print $$4; found_ua=1} } \
		{if ($$1 == "va" && found_va == 0) {print $$4; found_va=1} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'
