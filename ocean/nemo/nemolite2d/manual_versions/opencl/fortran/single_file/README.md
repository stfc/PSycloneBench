# OpenCL NEMOLite2D with Fortran driver layer #

This directory contains an OpenCL implementation of NEMOLite2D where
the driver layer is written in Fortran (using the clFortran interface
to OpenCL). The kernels must still be in OpenCL and are in the
../../../kernels/opencl directory.

## Compiling ##

The Makefile(s) pick up various settings from environment variables.
Please see the ../c/README.md file for details.

Note that although CC must be set to "g++" when compiling the OpenCL kernels,
leaving it set to this will cause the build of dl_timer to fail. It must
therefore be set to "gcc" once the kernels have been compiled.

## Running ##

Unlike the C version, this Fortran implementation reads the supplied
"namelist" file in order to get the details of the model run (size of
mesh, number of time steps etc.). In its current configuration the model
is intended mainly to be run in emulation mode in order to validate
the OpenCL implementation. Again, see the ../c/README.md file for
more details.

