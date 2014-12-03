# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

# Location of the code defining the GOcean API
API_DIR = ../api_v0.2
API_LIB = ${API_DIR}/gocean_api.a

EXECS = go2d

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = ${API_LIB}

all: $(EXECS)

# Normal targets
go2d: 
	${MAKE} MODULE_LIST="gocean2d_direct.o ${COMMON_MODULES}" go2d.exe

${API_LIB}: ${API_DIR}/*.?90
	${MAKE} -C ${API_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}" API_LIB="gocean_api.a"

gocean2d_direct.o: $(COMMON_MODULES)

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${API_DIR} -c $<

clean: 
	${MAKE} -C ${API_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90

allclean: clean
	rm -f *.exe fparser.log

docs:
	doxygen gocean2d.doxy.config
