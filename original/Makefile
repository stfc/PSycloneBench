include ../common/Makefile.include

# For timing module 
TIMING_DIR = ../common/

F90FLAGS += -I${TIMING_DIR}

LDFLAGS += ${TIMING_DIR}timing_mod.o ${TIMING_DIR}intel_timer_mod.o

EXECS = shallow_base

all: $(EXECS)

shallow_base.o: timing

timing:
	${MAKE} -C ${TIMING_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}"

%: %.o 
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	${MAKE} -C ${TIMING_DIR} clean
	rm -f *.o *.mod *.MOD $(EXECS)
