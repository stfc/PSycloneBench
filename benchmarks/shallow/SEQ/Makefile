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
INF_LIB = ${INF_DIR}/src/lib_fd.a

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

.PHONY: all shallow_gen shallow_base inf_lib timer_lib

all: $(EXECS)

# Targets involving the PSyclone code-generation framework
shallow_gen.f90: infrastructure_mod.f90 shallow_alg.f90
	psyclone -s ./ocl_transform.py -api gocean1.0 shallow_alg.f90 -oalg shallow_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as shallow_gen.f90
psy.f90: shallow_gen.f90

# The generated code depends on the generated Psy middle-layer
shallow_gen: inf_lib timer_lib
	${MAKE} MODULE_LIST="shallow_gen.o ${COMMON_MODULES} ${GENERATED_MODULES}" shallow_gen.exe

# Our module providing fake infrastructure kernels is in the INF directory
# but we need it to be in the same directory as all of the other kernel
# files.
infrastructure_mod.f90:
	ln -sf ${INF_INC}/infrastructure_mod.f90 .

# Normal targets
shallow_base: inf_lib timer_lib
	${MAKE} MODULE_LIST="shallow_base.o ${COMMON_MODULES}" shallow_base.exe

shallow_base.o: $(COMMON_MODULES)

shallow_gen.o: $(COMMON_MODULES) ${GENERATED_MODULES}

# Interdependencies between modules, alphabetical order

apply_bcs_mod.o: inf_lib
compute_cu_mod.o: inf_lib
infrastructure_mod.o: inf_lib
initial_conditions_mod.o: physical_params_mod.o model_mod.o
model_mod.o: inf_lib \
             physical_params_mod.o \
             shallow_io_mod.o \
             time_smooth_mod.o
physical_params_mod.o:  inf_lib
shallow_io_mod.o: inf_lib
time_smooth_mod.o: inf_lib
time_step_mod.o: inf_lib \
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

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

libclean:
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90 shallow_gen.f90
	rm -f gnu_opt_report.txt

allclean: clean libclean
	rm -f *.exe
	rm -rf *_wpl_dir
	rm -f infrastructure_mod.f90

docs:
	doxygen shallow.doxy.config
