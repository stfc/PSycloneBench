include ../psykal_vanilla/Makefile.include

# For timing module in PSyKAl version of code
TIMING_DIR = ../psykal_vanilla/

F90FLAGS += -I${TIMING_DIR}

LDFLAGS += ${TIMING_DIR}timing_mod.o ${TIMING_DIR}intel_timer_mod.o

EXECS = shallow_base

all: $(EXECS)

shallow_base.o: ${TIMING_DIR}timing_mod.o

%: %.o 
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD $(EXECS)
