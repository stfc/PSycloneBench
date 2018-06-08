# Manual OpenACC version of PSyKAl NEMOLite2D #

This directory contains the manual OpenACC implementation of the
PSyKAl form of NEMOLite2D. It does not use the 'common' NEMOLite2D code
base (in ../../common) because some of the routines have been modified
to support having data in a remote memory space.

## Compiling ##

You will need a compiler with OpenACC support e.g. PGI or Cray. (As of
June 2018, gfortran support for OpenACC is still in development - see
https://gcc.gnu.org/wiki/OpenACC.)

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/pgi_acc.sh for sample
settings for the PGI compiler suite.
If you are using Bash and the PGI compiler then:

    > . ../../../../../compiler_setup/pgi_acc.sh
    > make

should build the `nemolite2d.exe` binary.

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

Execution is then just:

    > ./nemolite2d.exe

If you wish to confirm that the GPU is being used (or see timings) then
set PGI_ACC_TIME=1 before running.

