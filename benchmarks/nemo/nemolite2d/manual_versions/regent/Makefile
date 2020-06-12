F90 = gfortran

%.o : %.f90
	$(F90) -fPIC -c -o $@ $<

libgocean2d_io_mod.so : kind_params_mod.o gocean2d_io_mod.o
	$(F90) kind_params_mod.o gocean2d_io_mod.o -shared -o $@

clean: 
	rm *.o *.mod *.so
