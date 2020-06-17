# OpenCL version of NEMOLite2D #

This directory contains the C implementation of the OpenCL version
of NEMOLite2D. It was originally developed to use the AMD dispatch
library with an NVIDIA GPU as the compute device. It has subsequently
been adapted to use the Intel OpenCL toolset to target an FPGA as
the compute device.

## Compiling ##

Although the code is general-purpose OpenCL, the Makefile is currently
configured to use the Altera (now Intel) Off-line Compiler. Therefore
you will either need to install this or to edit the Makefile to match your
configuration.

For FPGA *emulation* using the Intel (Altera) tools, INTELFPGA_PRO
must be substituted with the location in which you've installed the
Intel Quartus Pro toolkit and NALLATECH with the location of the Board
Support Package:

    # Settings needed to compile the FPGA OpenCL version of NEMOLite2D
    export OPENCL_LIBS="-L<INTELFPGA_PRO>/17.1/hld/board/nalla_pcie/linux64/lib -L<INTELFPGA_PRO>/17.1/hld/host/linux64/lib -Wl,--no-as-needed -lalteracl -lnalla_pcie_mmd -lelf"
    export OPENCL_INCLUDE="-I<INTELFPGA_PRO>/17.1/hld/host/include -I<NALLATECH>/common/inc"
    export F90=gfortran
    export CC=g++
    export CFLAGS="-g -fpermissive -Wall"

    export CL_CONTEXT_EMULATOR_DEVICE_ALTERA=1

Typing 'make nemolite2d_fpga' should then build the kernels and driver
code.

Of course, if you want to run on the FPGA for real then you need to
(go to bed while the code compiles and) unset CL_CONTEXT_EMULATOR_DEVICE_ALTERA.

The following configuration uses the Gnu compiler and AMD's dispatch library,
but the Makefile needs to be edited for this to work:

    export CC=gcc
    export CFLAGS="-g -O0 -Wall"
    export LDFLAGS=-lm
    export OPENCL_LIBS="-L/opt/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL"
    export OPENCL_INCLUDE=-I/opt/AMDAPPSDK-3.0/include

## Running ##

The benchmark problem size, number of time-steps and whether or not to
perform timing/profiling are also controlled via environment variables:

    export NEMOLITE2D_NSTEPS=10
    export NEMOLITE2D_NY=255
    export NEMOLITE2D_NX=255
    export NEMOLITE2D_PROFILING=1

OpenCL configuration is also performed by environment variables. For
running in emulation mode you must set:

    export CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA=1

In an OpenCL application, the kernels are loaded at run-time rather
than being linked into the executable. When using an FPGA, these
kernels must be compiled into a single image (.aocx file) ahead of
time. The location of this image is then supplied at runtime via
another environment variable:

    export NEMOLITE2D_SINGLE_IMAGE=${PWD}/nemolite2d_kernels.aocx

If profiling is enabled then the OpenCL profiling API is used to report
the amount of time spent in each kernel. This information is printed
to stdout at the end of a run.

## OpenCL Device Selection ##

Currently the code queries some of the properties of all of the OpenCL
devices it discovers and outputs this information to stdout.
e.g. on a system containing two NVIDIA GPUs it produces:

    Have 1 platforms.
    Platform 0 (id=36806128) is: NVIDIA CUDA
    Have 2 devices
    Device 0 is: Tesla K20c, type=4, version=OpenCL 1.2 CUDA
                 double precision supported
                 max work group size = 1024
                 max size of each work item dimension: 1024 1024
                 max compute units = 13
    Device 1 is: Quadro K600, type=4, version=OpenCL 1.2 CUDA
                 double precision supported
                 max work group size = 1024
                 max size of each work item dimension: 1024 1024
                 max compute units = 1

The code is currently set-up to choose device 0 as the one upon which
to execute the OpenCL kernels. This may need to be changed (in `nemolite.c`)
depending on your system.
