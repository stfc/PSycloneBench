# This Makefile expects the following env. vars to be set:
#
# F90       - How to invoke the Fortran compiler
# F90FLAGS  - Flags to pass to the Fortran compiler
# OMPFLAGS  - Flags for compiling with OpenMP
# AR        - Command to use when creating an archive (.a)

# Location of the dl_timer and infrastucture code
TIMER_DIR = ../../../shared/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ../../../shared/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a

EXECS = nemolite2d

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = ${INF_LIB}

all: $(EXECS)

# Normal targets
.PHONY: nemolite2d
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

timer_lib:
	make -C ${TIMER_DIR} sm_lib

inf_lib:
	make -C ${INF_DIR}

nemolite2d.o: $(COMMON_MODULES)

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) $(OMPFLAGS) -o $@ $(MODULE_LIST) $(TIMER_LIB) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

clean: 
	rm -f *.o *.mod *.MOD *~

allclean: clean
	rm -f *.exe
	rm -rf *_wpl_dir
	rm -f *.opt *.optrpt gnu_opt_report.txt
