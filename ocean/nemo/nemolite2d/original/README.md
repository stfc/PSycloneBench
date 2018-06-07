# Original NEMOLite2D #

This directory contains the original, single source-file version of
NEMOLite2D for reference. Note that it has been enhanced slightly and
thus makes use of checksum and I/O routines from both the dl_esm_inf
library and the `../common` directory.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite. e.g. in Bash you can set-up and then
build by doing:

    > . ../../../../compiler_setup/gnu.sh
    > make

This should create the `nemolite2d.exe` binary.

## Running ##

This directory contains various versions of the `namelist`
configuration file with the domain size as a suffix. The code is
distributed with a sym-link to the `namelist.128` file but you can of
course change this, e.g.:

    > ln -sf namelist.256 namelist

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing this `namelist` file.

Since the executable is serial, it may simply be executed on the command
line like so:

    > ./nemolite2d.exe

If doing performance benchmarking you may wish to ensure that the single
process is bound to a single processor core (see e.g. taskset, likwid or
Intel's KMP_AFFINITY).

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.
