# Makefile for the PSyKAl version of Shallow.

# The two targets that this Makefile supports
#  - shallow_base is the version of the code with manual invokes
#  - shallow_gen uses PSyClone to generate the invokes
EXECS = shallow_base shallow_gen

# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../shared

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a

# The modules that are common to both targets
MODULES = physical_params_mod.o \
          model_mod.o shallow_io_mod.o \
          initial_conditions_mod.o compute_cu_mod.o compute_cv_mod.o \
          compute_z_mod.o compute_h_mod.o \
          compute_unew_mod.o compute_vnew_mod.o compute_pnew_mod.o \
          apply_bcs_mod.o time_smooth_mod.o time_step_mod.o

# The modules only used when doing code generation
GENERATED_MODULES = infrastructure_mod.o psy.o

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${INF_LIB} 

all: $(EXECS)

# Targets involving the PSyclone code-generation framework
shallow_gen.f90: infrastructure_mod.f90 shallow_alg.f90
	psyclone -api gocean1.0 shallow_alg.f90 -oalg shallow_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as shallow_gen.f90
psy.f90: shallow_gen.f90

# The generated code depends on the generated Psy middle-layer
.PHONY: shallow_gen
shallow_gen:
	${MAKE} MODULE_LIST="shallow_gen.o ${COMMON_MODULES} ${GENERATED_MODULES}" shallow_gen.exe

# Our module providing fake infrastructure kernels is in the INF directory
# but we need it to be in the same directory as all of the other kernel
# files.
infrastructure_mod.f90:
	ln -sf ${INF_INC}/infrastructure_mod.f90 .

# Normal targets
.PHONY: shallow_base
shallow_base: inf_lib
	${MAKE} MODULE_LIST="shallow_base.o ${COMMON_MODULES}" shallow_base.exe

shallow_base.o: $(COMMON_MODULES)

shallow_gen.o: $(COMMON_MODULES) ${GENERATED_MODULES}

# Interdependencies between modules, alphabetical order

apply_bcs_mod.o: ${INF_LIB}
compute_cu_mod.o: ${INF_LIB}
infrastructure_mod.o: ${INF_LIB}
initial_conditions_mod.o: physical_params_mod.o model_mod.o
model_mod.o: ${INF_LIB} \
             physical_params_mod.o \
             shallow_io_mod.o \
             time_smooth_mod.o
physical_params_mod.o:  ${INF_LIB}
time_smooth_mod.o: ${INF_LIB}
time_step_mod.o: ${INF_LIB} \
                 apply_bcs_mod.o \
                 compute_cu_mod.o \
                 compute_cv_mod.o \
                 compute_pnew_mod.o \
                 compute_unew_mod.o \
                 compute_vnew_mod.o \
                 compute_z_mod.o \
                 compute_h_mod.o

%.exe: $(MODULE_LIST)
	$(F90) $(OMPFLAGS) -o $@ $(MODULE_LIST) $(TIMER_LIB) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

.PHONY: inf_lib timer_lib

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90 shallow_gen.f90

allclean: clean
	rm -f *.exe
	rm -rf *_wpl_dir
	rm -f infrastructure_mod.f90
	rm -f gnu_opt_report.txt

docs:
	doxygen shallow.doxy.config
