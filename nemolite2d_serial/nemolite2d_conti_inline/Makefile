# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

# Location of the code defining the GOcean API
API_DIR = ../../api_v1.0
API_LIB = ${API_DIR}/gocean_api.a
# Location of dl_timer
TIMER_DIR = ../../dl_timer
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a

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
COMMON_MODULES = ${TIMER_LIB} $(MODULES) ${API_LIB}

all: $(EXECS)

# Normal targets
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}" AR="${AR}" ARFLAGS="${ARFLAGS}" API_LIB="gocean_api.a"

${TIMER_LIB}:
	${MAKE} -C ${TIMER_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}" sm_lib

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

boundary_conditions_mod.o: physical_params_mod.o ${API_LIB} model_mod.o
continuity_mod.o: model_mod.o ${API_LIB}
gocean2d_io_mod.o: ${API_LIB}
model_mod.o: ${API_LIB} gocean2d_io_mod.o
momentum_mod.o: model_mod.o physical_params_mod.o ${API_DIR}/kind_params_mod.o
time_step_mod.o: ${API_LIB} momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o
time_update_mod.o: model_mod.o ${API_LIB}

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${API_DIR} -I${TIMER_DIR}/src -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${API_DIR} -I${TIMER_DIR}/src -c $<

clean: 
	${MAKE} -C ${API_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90

allclean: clean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
