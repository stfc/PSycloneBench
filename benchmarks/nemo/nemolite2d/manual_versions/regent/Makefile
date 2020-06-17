F90 = gfortran

%.o : %.f90
	$(F90) -fPIC -c -o $@ $<

all: libgocean2d_io_mod.so
regent: nemolite2d.exe

libgocean2d_io_mod.so : kind_params_mod.o gocean2d_io_mod.o
	$(F90) kind_params_mod.o gocean2d_io_mod.o -shared -o $@

nemolite2d.exe: 
	SAVEOBJ=1 $(REGENT_SCRIPT) algorithm.rg

clean: 
	rm *.o *.mod *.so
