# OpenCL version of the NEMOLite2D benchmark. #

The directory containing this file contains C and Fortran OpenCL
versions of the PSyKAl NEMOLite2D benchmark. It has been created by
converting the Fortran version into C and then adding the necessary
OpenCL infrastructure. The Fortran version uses clFortran
(github.com/cass-support/clfortran) for the driver layer.

## Dependencies ##

You will need Make, Fortran and C compilers and a working OpenCL
installation.  This code has been tested using AMD's OpenCL dispatch
library with NVIDIA's OpenCL library as a back-end.

## Building ##

The Makefile picks up the compiler to use, the location of the
OpenCL libraries and any flags from environment variables. 
Examples of the non-OpenCL-specific settings can be found in the
../../../../../compiler_setup directory.

e.g. to use the Gnu compiler and AMD's dispatch library:

export CC=gcc
export CFLAGS="-g -O0 -Wall"
export LDFLAGS=-lm
export OPENCL_LIBS="-L/opt/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL"
export OPENCL_INCLUDE=-I/opt/AMDAPPSDK-3.0/include

For FPGA *emulation* using the Intel (Altera) tools:

export CC=g++
export CFLAGS="-g -fpermissive"

export FPGA_BOARD=nalla_pcie
export OPENCL_LIBS="-L${INTELFPGAOCLSDKROOT}/board/${FPGA_BOARD}/linux64/lib -L${INTELFPGAOCLSDKROOT}/host/linux64/lib -Wl,--no-as-needed -lalteracl -l${FPGA_BOARD}_mmd -lelf"
export OPENCL_INCLUDE=-I${INTELFPGAOCLSDKROOT}/host/include

export CL_CONTEXT_EMULATOR_DEVICE_ALTERA=1

Typing 'make nemolite2d_fpga' should then build the kernels and driver
code.

Of course, if you want to run on the FPGA for real then you need to
(go to bed while the code compiles and) unset CL_CONTEXT_EMULATOR_DEVICE_ALTERA.

## Configuration and running ##

By default, the application runs a 128x128 domain for 2000 time-steps.
These settings may be configured at run-time through the following
environment variables; NEMOLITE2D_NX, NEMOLITE2D_NY and NEMOLITE2D_NSTEPS.
Note that NEMOLITE2D_N{X,Y} + 1 must be a multiple of 64 because of the
use of work-group sizes for the Momentum kernels.

Rudimentary OpenCL profiling may be enabled by setting NEMOLITE2D_PROFILING.
(There is a performance penalty associated with this.) The results of this
profiling are printed to stdout.

Host-side timing is performed on the CPU. The results of this are written
to `timing_report_serial.dat`.

The same source code is used to build both the OpenCL kernels and
those that are run on the CPU. At the end of a run, checksums are
computed for the various prognostic fields (`ua`, `va` and `ssha`) and the
results compared. The final values of these fields are also written to
file (e.g. `ua_cpu.dat` and `ua_ocl.dat`) to permit visualisation with
e.g. gnuplot.

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
