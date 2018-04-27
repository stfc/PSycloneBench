# Makefile.include includes the settings for the PGI compiler since
# this code uses OpenACC
include Makefile.include

# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../../../shared

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a

# Location of the code defining the GOcean API

EXECS = nemolite2d

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = ${INF_LIB}

all: $(EXECS)

# Normal targets
.PHONY: nemolite2d
nemolite2d: timer_lib inf_lib
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)

timer_lib:
	make -C ${TIMER_DIR} sm_lib

inf_lib:
	make -C ${INF_DIR}

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) $(OMPFLAGS) -o $@ $(MODULE_LIST) $(INF_LIB) $(TIMER_LIB) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

timer_lib:
	make -C ${TIMER_DIR} sm_lib

inf_lib:
	make -C ${INF_DIR}

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~

allclean: clean
	rm -f *.exe
	rm -rf *_wpl_dir

