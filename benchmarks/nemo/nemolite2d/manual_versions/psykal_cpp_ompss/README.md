# C++ PSyKal OmpSs version of NEMOLite2D #

This directory contains the OmpSs implementation of the PSyKAl form of
NEMOLite2D that has had its PSy layers manually implemented in C++.

The Fortran(NemoLite2d) to C++(PSy-layer) interoperability is done
with an `iso_c_binding` interface for the function `c_invoke_time_step`
in the `time_step_mod.f90`.
When compiling this function can be linked with object files defining
it in C or C++ (if the function is marked as `extern 'C'`).

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/ompss.sh for sample
settings for the OmpSs compiler suite. This version is untested without
the OmpSs compiler (though in principle it could run with other OpenMP 
compilers).
If you are using Bash and the OmpSs compiler then:

    > . ../../../../../compiler_setup/ompss.sh
    > make

The available target is:

nemolite2d_ompss - full NEMOLite2D with the C++ PSy layer (`time_step.cpp`),
                   with OmpSs task pragmas applied to the for loops in the 
                   PSy-later (`time_step_ompss.cpp`)

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

The OmpSs version can be configured with the OmpSs environment
variables to tune the number of threads and core affinities, or other
runtime parameters.

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.


