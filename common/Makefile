OBJECTS = intel_timer_mod.o timing_mod.o

all: ${OBJECTS}

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

timing_mod.o: intel_timer_mod.o

clean:
	rm -f *.o *.mod *.MOD

