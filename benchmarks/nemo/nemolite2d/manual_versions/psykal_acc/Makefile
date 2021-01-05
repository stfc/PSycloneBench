# Makefile for the manual, PSyKAl, OpenACC version of nemolite2d.
# Note, this does not yet use the 'common' routines in ../../common
# because the versions here have been modified slightly to
# support having data in a remote memory space.

# This Makefile picks up the compiler and flags etc. from environment
# variables. Example settings for the PGI compiler can be found in 
# ../../../../../compiler_setup/pgi_acc.sh

# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../../../shared

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/libdl_timer_omp.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/lib_fd.a

COMMON_DIR = ../../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a

# The targets that this Makefile supports
EXECS = nemolite2d

# The modules that are common to both targets
MODULES = gocean2d_io_mod.o \
          time_step_mod.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${INF_LIB}

.PHONY: all timer_lib inf_lib nemolite2d

all: $(EXECS)

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

${COMMON_LIB}:
	${MAKE} -C ${COMMON_DIR}

# Normal targets
nemolite2d: ${COMMON_LIB} timer_lib inf_lib
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

boundary_conditions_mod.o: inf_lib
gocean2d_io_mod.o: inf_lib
infrastructure_mod.o: inf_lib
model_mod.o: inf_lib gocean2d_io_mod.o
time_step_mod.o: inf_lib
time_update_mod.o: inf_lib

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) $(OMPFLAGS) -o $@ $(MODULE_LIST) ${COMMON_LIB} $(TIMER_LIB) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	${MAKE} -C ${TIMER_DIR} clean
	rm -f *.o *.mod *.MOD *~

allclean: clean
	rm -f *.exe
	rm -f *.lst *.cg *.opt *.optrpt *.s
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
