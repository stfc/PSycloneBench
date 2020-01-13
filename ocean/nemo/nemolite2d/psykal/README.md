# PSyKAl NEMOLite2D with Code Generation #

This directory contains the PSyKAl version of NEMOLite2D that uses
PSyclone to generate (and transform) the Parallel System (PSy) layer.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite. In Bash you can source that file
to set things up for Gnu:

    > . ../../../../compiler_setup/gnu.sh

Since the Makefile uses PSyclone you must have PSyclone installed on your
system (see http://psyclone.readthedocs.io/en/stable/getting_going.html).

The Makefile supports two targets:

1. nemolite2d_gen - vanilla (serial) version
2. nemolite2d_gen_omp - OpenMP version with static scheduling

but just executing `make` will build both and should result in the
creation of two executables - `nemolite2d_gen.exe` and
`nemolite2d_gen_omp.exe`. (The PSyclone-generated PSy code for these
two binaries can be found in `psy.f90` and `psy_omp.f90`, respectively.)

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

Since `nemolite2d_gen.exe` is serial, it can just be executed on
the command line, e.g.:

    > ./nemolite2d_gen.exe

For the OpenMP version, you will need to set OMP_NUM_THREADS to the
desired number of threads to use.

Note that if you are benchmarking performance then you probably need
to ensure that the threads are bound sensibly to processor cores (see
e.g. taskset, likwid or Intel's KMP_AFFINITY).

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.

