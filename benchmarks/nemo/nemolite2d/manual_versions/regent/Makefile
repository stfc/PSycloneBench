# Makefile for the Regent version of NemoLite2D
#
# This Makefile requires the REGENT_SCRIPT environment variable to 
# point to the regent.py script
# As the regent version requires libgfortran during execution, 
# gfortran is a requirement for this version.


F90 = gfortran

#-fPIC is required to enable the creation of the shared library.
%.o : %.f90
	$(F90) -O3 -fPIC -c -o $@ $<

all: libgocean2d_io_mod.so
regent: nemolite2d.exe

#Specific rule to create the shared library required for the Regent version.
libgocean2d_io_mod.so : kind_params_mod.o gocean2d_io_mod.o
	$(F90) -O3 kind_params_mod.o gocean2d_io_mod.o -shared -o $@

nemolite2d.exe: algorithm.rg bc_flather.rg calculate_sea_surface_height_t.rg calculate_velocity_fields.rg \
                checksum.rg initialise_grid_points.rg model_init.rg model_write.rg read_config.rg  \
                read_namelist.rg update_sea_surface_t.rg update_values.rg update_velocity_boundary.rg
	SAVEOBJ=1 $(REGENT_SCRIPT) algorithm.rg

clean: 
	rm *.o *.mod *.so
