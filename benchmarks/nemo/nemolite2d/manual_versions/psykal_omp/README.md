# Manual OpenMP version of PSyKAl NEMOLite2D #

This directory contains a version of the PSyKAl form of NEMOLite2D
that has had its PSy layer (in `time_step_mod.f90`) manually
parallelised using OpenMP.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
If you are using Bash and the Gnu compiler then:

    > . ../../../../../compiler_setup/gnu.sh
    > make

should build the `nemolite2d.exe` binary.

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

The number of threads to use must be set via the OMP_NUM_THREADS
environment variable. All of the parallelised loops use the `RUNTIME`
schedule option which means that the schedule to use is picked up
at run-time from the OMP_SCHEDULE environment variable. The simplest
option is simply to set this to "static":

    > export OMP_NUM_THREADS=4
    > export OMP_SCHEDULE="static"
    > ./nemolite2d.exe

Note that for performance benchmarking you will probably need to ensure
that your OpenMP threads are sensibly bound to processor cores (e.g.
using taskset, likwid or Intel's KMP_AFFINITY).

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.


