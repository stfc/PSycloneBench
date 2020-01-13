# Original Shallow Benchmark #

This is the 'original' version of Shallow, obtained from Graham Riley
at the University of Manchester at the beginning of the NERC GOcean
project (2014).

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
If you are using Bash and the Gnu compiler then:

    > . ../../../../compiler_setup/gnu.sh
    > make

will build the code.

## Running ##

If the compilation is successful then a binary named `shallow_base` will
be created.

Model parameters (size of domain, number of time-steps, whether or not
and how often to do output) may be configured by editing the
`namelist` file.

Execution is then just `./shallow_base`. Note that if you are
benchmarking performance then you probably need to ensure that the
process is bound to a single processor core (see e.g. taskset, likwid
or Intel's KMP_AFFINITY).

## Output ##

If output is enabled in the namelist file then various fields are dumped
to disk during the model run. The created files are ASCII and formatted
for use with `splot` in gnuplot.

