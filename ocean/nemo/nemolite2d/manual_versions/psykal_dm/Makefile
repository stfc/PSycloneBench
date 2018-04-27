# Makefile for manual MPI version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=mpif90
# export F90FLAGS="-O3"
# export AR=ar

# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../../../shared

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a

# The targets that this Makefile supports
EXECS = nemolite2d

# The modules that are common to both targets
MODULES = model_mod.o boundary_conditions_mod.o \
          continuity_mod.o initialisation_mod.o \
          momentum_mod.o physical_params_mod.o \
          time_update_mod.o gocean2d_io_mod.o \
          time_step_mod.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${INF_LIB}

timer_lib:
	make -C ${TIMER_DIR} dm_lib

inf_lib:
	make -C ${INF_DIR}

all: $(EXECS)

# Normal targets
.PHONY: nemolite2d
nemolite2d: timer_lib inf_lib
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

boundary_conditions_mod.o: physical_params_mod.o inf_lib model_mod.o
continuity_mod.o: model_mod.o inf_lib
gocean2d_io_mod.o: inf_lib
model_mod.o: inf_lib gocean2d_io_mod.o
momentum_mod.o: model_mod.o physical_params_mod.o inf_lib
time_step_mod.o: inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o
time_update_mod.o: model_mod.o inf_lib

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) $(TIMER_LIB) $(LDFLAGS) ${OMPFLAGS}

%.o: %.f90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${INF_INC} -I${TIMER_INC} -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~

allclean: clean
	make -C ${TIMER_DIR} allclean
	make -C ${INF_DIR} distclean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
