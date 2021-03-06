# Makefile for manual MPI version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler and flags from
# environment variables.

# If F90 already has a custom MPI compiler command, it will use it,
# otherwise it will form a mpich or openmpi -family command using the
# compiler defined in F90.
ifndef F90
$(error The F90 env variable must be defined)
else
ifeq ($(findstring mpi,$(F90)),)  # If a 'mpi' string is not found in F90, the
	                              # Makefile will try to construct it itself.
# First we check that if the mpifort command is available in the environment.
ifneq ($(shell mpifort --version > /dev/null 2> /dev/null; echo $$?),0)
$(error This Makefile requires an MPI Fortran wrapper in the environment.)
else
# If it is available, we figure out which type of MPI implementation it is
# by testing the --showme flag (which is only available in OpenMPI).
ifeq ($(shell mpifort --showme > /dev/null 2>/dev/null; echo $$?),0)
# Construct F90 for OpenMPI-family MPI implementations.
F90 := OMPI_MPIF90=$(F90) mpifort
else
# Construct F90 for mpich-family MPI implementations.
F90 := mpifort -fc=${F90}
endif
endif
endif
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
	$(F90) -o $@ $(MODULE_LIST) $(TIMER_LIB_DM) $(LDFLAGS) ${OMPFLAGS}

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
