# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

# Location of the psyclone source files (for doing the code
# generation)
PSYCLONE_DIR = ../../PSyclone/src

# Location of the code defining the GOcean API
API_DIR = ../api_v0.2
API_LIB = ${API_DIR}/gocean_api.a

# The two targets that this Makefile supports
#  - shallow_base is the version of the code with manual invokes
#  - shallow_gen uses PSyClone to generate the invokes
EXECS = gocean2d

# The modules that are common to both targets
MODULES = model_mod.o momentum_mod.o physical_params_mod.o

GENERATED_MODULES = psy.o

COMMON_MODULES = ${API_LIB} $(MODULES)

all: $(EXECS)

# Targets involving the code-generation framework
shallow_gen.f90: shallow_gocean.f90
	python ${PSYCLONE_DIR}/generator.py -api gocean shallow_gocean.f90 -oalg shallow_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as shallow_gen.f90
psy.f90: shallow_gen.f90

# The generated code depends on the generated Psy middle-layer
shallow_gen:
	@echo "Code generation not yet working for version 0.2 of the GOcean API!"
	#${MAKE} MODULE_LIST="${COMMON_MODULES} ${GENERATED_MODULES}" shallow_gen.exe

# Normal targets
gocean2d: 
	${MAKE} MODULE_LIST="gocean2d_direct.o ${COMMON_MODULES}" go2d.exe

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}" API_LIB="gocean_api.a"

gocean2d_direct.o: $(COMMON_MODULES)

shallow_gen.o: $(COMMON_MODULES) ${GENERATED_MODULES}

# Interdependencies between modules, alphabetical order

momentum_mod.o: model_mod.o physical_params_mod.o ${API_DIR}/kind_params_mod.o

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

clean: 
	${MAKE} -C ${API_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90

allclean: clean
	rm -f *.exe fparser.log

docs:
	doxygen shallow.doxy.config
