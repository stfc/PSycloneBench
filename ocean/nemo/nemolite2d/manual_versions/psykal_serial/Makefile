# This Makefile expects the following env. vars to be set:
#
# F90       - How to invoke the Fortran compiler
# F90FLAGS  - Flags to pass to the Fortran compiler
# OMPFLAGS  - Flags for compiling with OpenMP
# AR        - Command to use when creating an archive (.a)

# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../../../shared

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a
COMMON_DIR = ../../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a

LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}

# The kernels used by this application and their location
KERNEL_DIR = ../kernels/fortran
KERNELS = boundary_conditions_mod.o \
          continuity_mod.o \
          momentum_mod.o \
          time_update_mod.o \
          infrastructure_mod.o

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(KERNELS) ${COMMON_LIB} ${INF_LIB}

# The targets that this Makefile supports

all: nemolite2d nemolite2d_align nemolite2d_cont_only nemolite2d_nodiv \
nemolite2d_preload nemolite2d_alignall nemolite2d_nopeel nemolite2d_8arrays

# Normal targets. Need to link with OpenMP in order to use OMP timer
nemolite2d: inf_lib ${KERNELS} time_step_mod.o nemolite2d.o
	${MAKE} MODULE_LIST="nemolite2d.o time_step_mod.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d.exe

# Version of nemolite2d with the time-stepping routine hacked to
# encourage the Intel compiler to vectorise efficiently
nemolite2d_align: ${COMMON_MODULES} time_step_mod_align.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_align.o nemolite2d.o ${COMMON_MODULES} ${LIBS}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_align.exe

# Time-stepping routine hacked to only call the Continuity kernel
nemolite2d_cont_only: ${COMMON_MODULES} time_step_mod_cont_only.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_cont_only.o nemolite2d.o ${COMMON_MODULES} ${LIBS}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_cont_only.exe

# With division ops removed from Momentum kernels
nemolite2d_nodiv: ${COMMON_MODULES} time_step_mod_align_nodiv.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_nodiv
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_align_nodiv.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_nodiv.exe

nemolite2d_preload: ${COMMON_MODULES} time_step_mod_align_preload.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_align_preload.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_preload.exe

nemolite2d_alignall: ${COMMON_MODULES} time_step_mod_alignall_preload.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod*.o time_step_mod*.mod
	${MAKE} MODULE_LIST="time_step_mod_alignall_preload.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_alignall.exe

nemolite2d_nopeel: ${COMMON_MODULES} time_step_mod_alignall_preload_nopeel.o \
nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod*.o time_step_mod*.mod
	${MAKE} MODULE_LIST="time_step_mod_alignall_preload_nopeel.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_nopeel.exe

nemolite2d_8arrays: ${COMMON_MODULES} time_step_mod_alignall_preload_nopeel_8arrays.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod*.o time_step_mod*.mod
	${MAKE} MODULE_LIST="time_step_mod_alignall_preload_nopeel_8arrays.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_8arrays.exe

# Interdependencies between modules, alphabetical order

time_step_mod.o: inf_lib ${KERNELS}

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) ${TIMER_LIB} $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<
%.o: %.F90
	$(F90) $(F90FLAGS) -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

# If we need a .f90 file that doesn't exist then it must be a kernel.
# Create a link to the required file...
%.f90:
	ln -sf ../../kernels/fortran/$@ .

.PHONY: inf_lib timer_lib

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

${COMMON_LIB}:
	${MAKE} -C ${COMMON_DIR}

clean: 
	rm -f *.o *.mod *.MOD *~

libclean:
	${MAKE} -C ${INF_DIR} clean

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log
	rm -f *.lst *.cg *.opt *.optrpt gnu_opt_report.txt *.s
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
