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
INF_LIB = ${INF_DIR}/src/lib_fd.a

COMMON_DIR = ../../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a

# The targets that this Makefile supports
EXECS = nemolite2d
MODULES = time_step_mod.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES)

.PHONY: all nemolite2d timer_lib inf_lib

all: $(EXECS)

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

${COMMON_LIB}:
	${MAKE} -C ${COMMON_DIR}

# Normal targets
nemolite2d: ${COMMON_LIB} timer_lib inf_lib time_step.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step.o" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

time_step_mod.o: ${COMMON_LIB} inf_lib

# Generic rules

# It needs the -stdc++ to link also the std lib functions called from C++ code
%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) $(LDFLAGS) ${OMPFLAGS} -lstdc++

time_step.o: time_step.cpp
	$(CXX) $(CFLAGS) -c $<

%.o: %.f90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.dat
	rm -f gnu_opt_report.txt

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

docs:
	doxygen gocean2d.doxy.config
