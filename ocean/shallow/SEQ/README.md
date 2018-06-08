# PSyKAl version of Shallow Benchmark #

This directory contains the Shallow benchmark obtained by re-structuring
to conform to the PSyKAl (Parallel System, Kernel, Algorithm) separation
of concerns.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
If you are using Bash and the Gnu compiler then:

    > . ../../../compiler_setup/gnu.sh
    > make shallow_base

will build the manual version of the code and produce `shallow_base.exe`.

The Makefile also supports a second target, 'shallow_gen', which uses
PSyclone to generate the Parallel-System (PSy) layer. For this to work
your system must have PSyclone installed
(see http://psyclone.readthedocs.io/en/stable/getting_going.html).
Once you have PSyclone, doing:

    > make shallow_gen

should create `shallow_gen.exe`.

## Running ##

Model parameters (size of domain, number of time-steps, whether or not
and how often to do output) may be configured by editing the
`namelist` file.

Execution is then just `./shallow_base.exe` or
`./shallow_gen.exe`. Note that if you are benchmarking performance
then you probably need to ensure that the process is bound to a single
processor core (see e.g. taskset, likwid or Intel's KMP_AFFINITY).

## Output ##

If output is enabled in the namelist file then various fields are dumped
to disk during the model run. The created files are ASCII and formatted
for use with `splot` in gnuplot.
