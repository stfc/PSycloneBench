# Makefile for the manual OpenACC version of nemolite2d.
# Makefile.include contains the settings for the PGI compiler
include ../manual_acc/Makefile.include

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
MODULES = model_mod.o \
          initialisation_mod.o \
          physical_params_mod.o \
          gocean2d_io_mod.o \
          time_step_mod.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${INF_LIB}

all: $(EXECS)

timer_lib:
	make -C ${TIMER_DIR} sm_lib

inf_lib:
	make -C ${INF_DIR}

# Normal targets
.PHONY: nemolite2d
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

boundary_conditions_mod.o: physical_params_mod.o inf_lib model_mod.o
gocean2d_io_mod.o: inf_lib
infrastructure_mod.o: inf_lib
model_mod.o: inf_lib gocean2d_io_mod.o
time_step_mod.o: inf_lib physical_params_mod.o
# boundary_conditions_mod.o
time_update_mod.o: model_mod.o inf_lib

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) $(OMPFLAGS) -o $@ $(MODULE_LIST) $(TIMER_LIB) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~

allclean: clean
	rm -f *.exe
	rm -f *.lst *.cg *.opt *.optrpt *.s
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
