include ../Makefile.include

API_DIR = ../../../api_v0.2
API_LIB = ${API_DIR}/gocean_api.a

# For timing module in new version of code
F90FLAGS += -I${API_DIR}

LDFLAGS += ${API_LIB} -lnetcdf -lnetcdff

EXECS = shallow_base

all: $(EXECS)

shallow_base.o: ${API_LIB}

%: %.o 
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD $(EXECS)
