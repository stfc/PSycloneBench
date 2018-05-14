# Makefile for manual OpenMP version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=gfortran
# export F90FLAGS="-O3"
# export OMPFLAGS=" -fopenmp"
# export AR=ar

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

# The targets that this Makefile supports
EXECS = nemolite2d

KERNELS = boundary_conditions_mod.o \
          continuity_mod.o \
          momentum_mod.o \
          time_update_mod.o \
          infrastructure_mod.o

MODULES = time_step_mod.o ${KERNELS}

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${INF_LIB}

.PHONY: all nemolite2d timer_lib inf_lib

all: $(EXECS)

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

${COMMON_LIB}:
	${MAKE} -C ${COMMON_DIR}

# Normal targets
nemolite2d: ${COMMON_LIB} timer_lib inf_lib kernels
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d.o: kernels $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

time_step_mod.o: ${COMMON_LIB} inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB} $(TIMER_LIB) $(LDFLAGS) ${OMPFLAGS}

%.o: %.f90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

kernels:
	ln -sf ../../kernels/fortran/*.?90 .

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
