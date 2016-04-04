include ../Makefile.include

API_VERSION=1.0
API_DIR = ../../../api_v${API_VERSION}
API_LIB = ${API_DIR}/gocean_api.a

# For timing module in new version of code
F90FLAGS += -I${TIMER_INC}

# For IO module
F90FLAGS += -I../

LDFLAGS += ${API_LIB} -lnetcdf -lnetcdff

EXECS = shallow_base

all: $(EXECS)

shallow_base.o: ${API_LIB} ../shallow_io_mod.o

../shallow_io_mod.o:
	${MAKE} -C ../

%: %.o 
	$(F90) $(OMPFLAGS) $(F90FLAGS) -o $@ $^ ../shallow_io_mod.o $(TIMER_LIB) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD $(EXECS)
