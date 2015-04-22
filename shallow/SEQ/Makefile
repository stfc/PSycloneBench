# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file in the arch dir.
include Makefile.include

# Location of the psyclone source files (for doing the code
# generation)
PSYCLONE_DIR = ../../../PSyclone/src

# Location of the code defining the GOcean API
API_DIR = ../../api_v0.2
API_LIB = ${API_DIR}/gocean_api.a

# The two targets that this Makefile supports
#  - shallow_base is the version of the code with manual invokes
#  - shallow_gen uses PSyClone to generate the invokes
EXECS = shallow_base shallow_gen

# The modules that are common to both targets
MODULES = physical_params_mod.o \
          model_mod.o shallow_io_mod.o \
          initial_conditions_mod.o compute_cu_mod.o compute_cv_mod.o \
          compute_z_mod.o compute_h_mod.o \
          compute_unew_mod.o compute_vnew_mod.o compute_pnew_mod.o \
          apply_bcs_mod.o time_smooth_mod.o time_step_mod.o

GENERATED_MODULES = psy.o shallow_gen.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${API_LIB} 

all: $(EXECS)

# Targets involving the code-generation framework
shallow_gen.f90: shallow_gocean.f90
	python ${PSYCLONE_DIR}/generator.py -api gocean1.0 shallow_gocean.f90 -oalg shallow_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as shallow_gen.f90
psy.f90: shallow_gen.f90

# The generated code depends on the generated Psy middle-layer
shallow_gen:
	${MAKE} MODULE_LIST="${COMMON_MODULES} ${GENERATED_MODULES}" shallow_gen.exe

# Normal targets
shallow_base: ${API_LIB}
	${MAKE} MODULE_LIST="shallow_base.o ${COMMON_MODULES}" shallow_base.exe

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}" API_LIB="gocean_api.a"

shallow_base.o: $(COMMON_MODULES)

shallow_gen.o: $(COMMON_MODULES) ${GENERATED_MODULES}

# Interdependencies between modules, alphabetical order

compute_cu_mod.o: ${API_LIB}

initial_conditions_mod.o: physical_params_mod.o model_mod.o

time_smooth_mod.o: ${API_LIB}

time_step_mod.o: ${API_LIB} \
                 apply_bcs_mod.o \
                 compute_cu_mod.o \
                 compute_cv_mod.o \
                 compute_pnew_mod.o \
                 compute_unew_mod.o \
                 compute_vnew_mod.o \
                 compute_z_mod.o \
                 compute_h_mod.o

model_mod.o: ${API_LIB} \
             physical_params_mod.o \
             shallow_io_mod.o \
             time_smooth_mod.o

apply_bcs_mod.o: ${API_LIB}

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

clean: 
	${MAKE} -C ${API_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90 shallow_gen.f90

allclean: clean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir

docs:
	doxygen shallow.doxy.config
