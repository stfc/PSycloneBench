# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file in the arch dir.
include Makefile.include

# Location of the psyclone source files (for doing the code
# generation)
PSYCLONE_DIR = ../../../PSyclone/src

# Location of the code defining the GOcean API
API_VERSION=1.0
API_DIR = ../../api_v${API_VERSION}
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

# The modules only used when doing code generation
GENERATED_MODULES = infrastructure_mod.o psy.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${API_LIB} 

all: $(EXECS)

# Targets involving the code-generation framework
shallow_gen.f90: infrastructure_mod.f90 shallow_alg.f90
	python ${PSYCLONE_DIR}/generator.py -api gocean1.0 shallow_alg.f90 -oalg shallow_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as shallow_gen.f90
psy.f90: shallow_gen.f90

# The generated code depends on the generated Psy middle-layer
shallow_gen:
	${MAKE} MODULE_LIST="shallow_gen.o ${COMMON_MODULES} ${GENERATED_MODULES}" shallow_gen.exe

# Our module providing fake infrastructure kernels is in the API directory
# but we need it to be in the same directory as all of the other kernel
# files.
infrastructure_mod.f90:
	ln -sf ${API_DIR}/infrastructure_mod.f90 .

# Normal targets
shallow_base: ${API_LIB}
	${MAKE} MODULE_LIST="shallow_base.o ${COMMON_MODULES}" shallow_base.exe

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} AR="${AR}" F90="${F90}" F90FLAGS="${F90FLAGS}" API_LIB="gocean_api.a"

shallow_base.o: $(COMMON_MODULES)

shallow_gen.o: $(COMMON_MODULES) ${GENERATED_MODULES}

# Interdependencies between modules, alphabetical order

apply_bcs_mod.o: ${API_LIB}
compute_cu_mod.o: ${API_LIB}
infrastructure_mod.o: ${API_LIB}
initial_conditions_mod.o: physical_params_mod.o model_mod.o
model_mod.o: ${API_LIB} \
             physical_params_mod.o \
             shallow_io_mod.o \
             time_smooth_mod.o
physical_params_mod.o:  ${API_LIB}
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
