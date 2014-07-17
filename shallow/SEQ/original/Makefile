include ../Makefile.include

# For timing module in new version of code
F90FLAGS += -I../

EXECS = shallow_base

OBJECTS = tiling_mod.o ../timing_mod.o ../intel_timer_mod.o

all: $(EXECS)

shallow_base.o: ../timing_mod.o tiling_mod.o

%: %.o 
	$(F90) $(F90FLAGS) $(OBJECTS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD $(EXECS)
