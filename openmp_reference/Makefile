include ../common/Makefile.include

# For timing module
TIMING_DIR = ../common/

F90FLAGS += -I${TIMING_DIR} ${OMPFLAGS}

LDFLAGS += ${TIMING_DIR}timing_mod.o ${TIMING_DIR}intel_timer_mod.o

EXECS = shallow_omp

all: $(EXECS)

shallow_omp.o: timing

%: %.o
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

timing:
	${MAKE} -C ${TIMING_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}"

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

clean: 
	${MAKE} -C ${TIMING_DIR} clean
	rm -f *.o *.mod *.MOD $(EXECS)
