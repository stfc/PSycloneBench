# Manual Kokkos version of PSyKAl NEMOLite2D #

This directory contains a version of the PSyKAl form of NEMOLite2D
that has had its PSy layer (in `time_step_kokkos.cpp`) manually
parallelised using Kokkos. The Fortran(NemoLite2d) to C++(PSy-layer)
is done with the interface introduced in the `psykal_cpp` version.

Kokkos is a Performance Portability EcoSystem for C++,
a good introduction to its concepts and APIs can be found in their
(wiki)[https://github.com/kokkos/kokkos/wiki].


## Implementation limitations ##

Kokkos has several abstractions to control where (execution spaces) and
how (execution policies and execution patterns) computation is executed
and how memory is laid out (Views).

This initial implementation is limited to a single execution space (OpenMP)
and is not yet using Views (essential for performance portability).

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
Additionally, it needs the path to the Kokkos installation, this can be
set inside the makefile or with the `KOKKOS_PATH` environment variable.
If you are using Bash and the Gnu compiler then:

    > . ../../../../../compiler_setup/gnu.sh
    > export KOKKOS_PATH=<path to kokkos>
    > make

should build the `nemolite2d_kokkos.exe` binary.

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

Kokkos supports multiple execution spaces (e.g Serial, Threads, OpenMP,
Cuda, HPX) but the current implementation is limited to only OpenMP.

Therefore, we can configure the parallel environment using the OpenMP
environment variables. For example, the number of threads to use can be
set via the OMP_NUM_THREADS environment variable.

    > export OMP_NUM_THREADS=4
    > ./nemolite2d_kokkos.exe

This will change in follow up implementations.

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.
