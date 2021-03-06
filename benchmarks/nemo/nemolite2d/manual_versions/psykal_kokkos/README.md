# Manual Kokkos version of PSyKAl NEMOLite2D #

This directory contains a version of the PSyKAl form of NEMOLite2D
that has had its PSy layer manually parallelised using Kokkos.
The Fortran(NemoLite2d) to C++(PSy-layer) is done with the interface
introduced in the `psykal_cpp` version.

Kokkos is a Performance Portability EcoSystem for C++,
a good introduction to its concepts and APIs can be found in their
(wiki)[https://github.com/kokkos/kokkos/wiki].


## Compiling and Implementations ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite. Additionally, it needs the path to
the Kokkos source code, this can be set inside the makefile or with the
`KOKKOS_PATH` environment variable. Note that the linking step is done
with the C++ compiler, but as it links also Fortran files, the libgfortran
library needs to be in the library path. If you are using Bash and the Gnu
compiler then:

    > source ../../../../../compiler_setup/gnu.sh
    > export KOKKOS_PATH=<path to kokkos>
    > make

will compile the default `nemolite2d_views_kokkos_OpenMP.exe` binary.

Kokkos supports multiple execution spaces (e.g Serial, Threads, OpenMP,
Cuda, HPX), these can be set for each Makefile target with the
KOKKOS_DEVICES environment variable (examples below). If none is provided,
it defaults to the OpenMP backend.

Kokkos also has several abstractions to control where (execution spaces)
and how (execution policies and execution patterns) computation is executed
and how memory is laid out (Views). This make it possible to create very
different Kokkos implementations, at the moment there are 2 versions:

- Rawpointers version: This version uses the Kokkos parallel dispatch
(execution policies and execution pattern) but on top of a raw pointer
arrays (given by the dl_esm_inf infrastructure). This
version is memory space efficient as it doesn't need to do extra copies
between the Fortran and Kokkos parts of the code, but it is not able to
run on the GPU or change the data layout. This version is available in
`time_step_rawpointers_kokkos.cpp` and can be build with the command:

    > make nemolite2d_rawpointers_kokkos

- View Containers: This version uses the Kokkos View containers in addition
to the Kokkos parallel dispatch. This allows Kokkos to control the data layout,
the padding, and the synchonization between host and device (GPU execution)
but it requires to keep two copies of the simulation data.
This version is available in `time_step_views_kokkos.cpp` and can be built
with an OpenMP or a Cuda backend by setting the KOKKOS_DEVICES environment
variable. Note that the Cuda back-end requires that the `nvcc` compiler is
installed on the system and available in PATH. See below examples of how to
compile the Kokkos View version for different devices:

    > make nemolite2d_views_kokkos KOKKOS_DEVICES=OpenMP

    > make nemolite2d_views_kokkos KOKKOS_DEVICES=Cuda

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

Then the NemoLite2D simulation is started by executing the chosen binary
(with the desired implementation and execution device) previously produced
by the build system:

    > ./nemolite2d_<version>_<device>.exe

In the case of OpenMP we can configure the parallel environment using the
standard OpenMP environment variables. For example, the number of threads
to use can be set via the OMP_NUM_THREADS environment variable.

    > export OMP_NUM_THREADS=4
    > ./nemolite2d_views_kokkos_OpenMP.exe


## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.
