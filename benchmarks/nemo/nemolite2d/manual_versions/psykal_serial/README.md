# Manual, serial, PSyKal versions of NEMOLite2D #

This directory contains various versions of the manual PSy layer
intended to be used to explore the performance of the benchmark.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
If you are using Bash and the Gnu compiler then:

    > . ../../../../../compiler_setup/gnu.sh
    > make

should build all of the various versions:

nemolite2d.exe - full NEMOLite2D with all kernels in-lined into the PSy
                 layer (time_step_mod.f90).
nemolite2d_align.exe - as above but with modifications to ensure the Intel
                       compiler knows arrays are aligned in the Continuity
		       kernel.
nemolite2d_cont_only.exe - PSy layer hacked so only contains Continuity kernel.
nemolite2d_nodiv.exe - as nemolite2d_align but with the Continuity kernel
                       hacked to remove division operations.
nemolite2d_preload.exe - as above but with the invoke of the Continuity kernel
                         modified to ensure all arrays are loaded into cache
			 before timing commences.
nemolite2d_alignall.exe - as above but with the Continuity kernel further
                          modified so that all array accesses are aligned.
nemolite2d_nopeel.exe - as above but with the explicit peel loop removed from
		        the invoke of the Continuity kernel.
nemolite2d_8arrays.exe - as above but with the Continuity kernel hacked
                         further so that only 8 arrays are accessed.

## Running ##

Since the various executables are all serial, they can just be executed on
the command line, e.g.:

    > ./nemolite2d.exe

Note that if you are benchmarking performance then you probably need
to ensure that the process is bound to a single processor core (see
e.g. taskset, likwid or Intel's KMP_AFFINITY).

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.


