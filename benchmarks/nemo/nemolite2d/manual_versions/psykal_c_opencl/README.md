# PSyKAl C OpenCL version of NEMOLite2D #

This directory contains an OpenCL implementation of NEMOLite2D that
uses the PSyKAl form with a Fortran to C first PSy-layer and a pure
C OpenCL second PSy-layer.

## Compiling ##

This implementation supports both OpenCL models: JIT (e.g. Intel OpenAPI
CPUs, NVidia OpenCL, AMD OpenCL, ... ) and ahead-of-time compilation (e.g.
Xilinx OpenCL, Intel OneAPI FPGAs, ...).

The host program is the same for both, and it is compiled with the default
Makefile target:

    make

this will generate a `nemlite2d.exe` file.

There are additional Makefile targets to perform the ahead-of-time/offline
compilation of supported compiler tool-chains. These are:

 - A target for the Xilinx OpenCl `allkernels_xilinx.xclbin` binary.

    make xilinxoffline

 - A target for the Intel Altera `allkernels_altera.aocx` binary.

    make alteraoffline

See the section below for instructions to provide the generated OpenCL
kernels binaries to the host OpenCL application.

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps [nitend],
whether or not and how often to do output [irecord]) may be configured by editing
the namelist file.

Additionally there are some environment parameters to control the OpenCL
execution, these are:

- `NEMOLITE2D_SINGLE_IMAGE`: When using ahead-of-time compied kernels this can
be provided to NemoLite2d as a single binary image in this environment variable.

- `OPENCL_PLATFORM`: Optional integer 0-base index number to specify in which OpenCL
platform the kernels should be launched on. If not specified it launches them in
platform 0. Use `clinfo` command to list the available platforms in the system.

- `NEMOLITE2D_PROFILING`: If this environment variable exist, it will set up the
OpenCL queue with the `CL_QUEUE_PROFILING_ENABLE` property.

Note that for running in each of the OpenCL devices the system should be properly
configured with the OpenCL ICDs and JIT compilers necessary. For instance, in the
Hartree Centre Wheatley workstation, the following set of commands can be executed
to run the code in all available platforms:

    module load NVidiaCuda
    ./nemolite2d.exe # Will run on the GPU
    module unload NVidiaCuda

    module load IntelOneAPI
    ./nemolite2d.exe # Will run on the FPGA software emulator
    OPENCL_PLATFORM=1 ./nemolite2d.exe # Will run on the multi-core Xeon CPU
    module unload IntelOneAPI

    module load XilinxVitis
    NEMOLITE2D_SINGLE_IMAGE=allkernels_xilinx.xclbin ./nemolite2d.exe # Will run on the FPGA
    module unload XilinxVitis

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.

