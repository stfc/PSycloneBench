# Makefile for the 'original', single-file version of Shallow.
# You will need to set which compiler to use etc. via environment
# variables. See e.g. ../../../../compiler_setup/gfortran.sh

SHARED_DIR=../../../../shared
TIMER_DIR=${SHARED_DIR}/dl_timer
TIMER_INC=${TIMER_DIR}/src
TIMER_LIB=${TIMER_DIR}/dl_timer_lib.a

# For timing module in new version of code
F90FLAGS += -I${TIMER_INC}

# For IO module
F90FLAGS += -I../

EXECS = shallow_base

.PHONY: all timer_lib libclean

all: $(EXECS)

shallow_base.o: timer_lib ../shallow_io_mod.o

../shallow_io_mod.o:
	${MAKE} -C ../

%: %.o 
	$(F90) $(OMPFLAGS) $(F90FLAGS) -o $@ $^ ../shallow_io_mod.o $(TIMER_LIB) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

libclean:
	${MAKE} -C ${TIMER_DIR} allclean

clean: 
	rm -f *.o *.mod *.MOD gnu_opt_report.txt

allclean: clean libclean
	rm -f $(EXECS)
