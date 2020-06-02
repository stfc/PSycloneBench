# C++ PSyKal versions of NEMOLite2D #

This directory contains various version of the PSyKAl form of
NEMOLite2D that has had its PSy layers manually implemented in C++.

The Fortran(NemoLite2d) to C++(PSy-layer) interoperability is done
with a `iso_c_binding` interface of the function `c_invoke_time_step`
in the `time_step_mod.f90`.
When compiling this function can be linked with object files defining
it in C or C++ (if the function is marked as `extern 'C'`).

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
If you are using Bash and the Gnu compiler then:

    > . ../../../../../compiler_setup/gnu.sh
    > make

The available targets are:

nemolite2d - full NEMOLite2D with the C++ PSy layer (`time_step.cpp`).
nemolite2d_omp - as above but with work-sharing OpenMP pragmas in for each
                 loop in the PSy-layer (`time_step_omp.cpp`). 
nemolite2d_ompopt - as above but with some optimizations to achieve better
                    performance (`time_step_ompopt.cpp`).

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

The versions with OpenMP can be configured with the OpenMP environment
variables to tune the number of threads and core affinities.

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.


