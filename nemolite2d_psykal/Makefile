# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

# Location of the psyclone source files (for doing the code
# generation)
PSYCLONE_DIR = ../../PSyclone/src

# Location of the code defining the GOcean API
API_DIR = ../api_v0.2
API_LIB = ${API_DIR}/gocean_api.a

# The targets that this Makefile supports
EXECS = nemolite2d nemolite2d_gen

# The modules that are common to both targets
MODULES = model_mod.o boundary_conditions_mod.o \
          continuity_mod.o initialisation_mod.o \
          momentum_mod.o physical_params_mod.o \
          time_update_mod.o gocean2d_io_mod.o \
          time_step_mod.o

GENERATED_MODULES = psy.o nemolite2d_gen.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${API_LIB}

all: $(EXECS)

# Targets involving the code-generation framework
nemolite2d_gen.f90: nemolite2d_alg.f90
	python ${PSYCLONE_DIR}/generator.py -api gocean1.0 nemolite2d_alg.f90 -oalg nemolite2d_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as nemolite2d_gen.f90
psy.f90: nemolite2d_gen.f90

# The generated code depends on the generated Psy middle-layer
nemolite2d_gen:
	${MAKE} MODULE_LIST="${GENERATED_MODULES} ${COMMON_MODULES}" nemolite2d_gen.exe

# Normal targets
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}" AR="${AR}" ARFLAGS="${ARFLAGS}" API_LIB="gocean_api.a"

nemolite2d.o: $(COMMON_MODULES)

nemolite2d_gen.o: $(COMMON_MODULES) psy.o

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
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

clean: 
	${MAKE} -C ${API_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90 nemolite2d_gen.f90

allclean: clean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
