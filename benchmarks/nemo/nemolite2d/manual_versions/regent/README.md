# Regent version of NEMOLite2D #

This implementation of NEMOLite2D uses Regent (http://regent-lang.org/) to handle the dataflow parallemlism.

Requirements:
LLVM 3.5-3.8 or 6
Python 2.7 or 3.X

Options:
GASNet (for distributed memory, yet to be tested)
CUDA 5.0 or later (for GPUs, yet to be implemented/tested)

## Building ##
The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
If you are using Bash and the Gnu compiler then:

    > . ../../../../../compiler_setup/gnu.sh
    > make

should build the required library for the Regent example (libgocean2d_io_mod.so)

To build a binary version of the Regent code (optional), run

    > export REGENT_SCRIPT=/patch/to/regent.py
    > make regent

## Running ##

If you built the binary version of the Regent code, you can just run that executable
provided the libgocean2d_io_mod.so and libregent.so libraries are visible on $LD_LIBRARY_PATH.

You can also run the code with the just-in-time regent compiler, using

    > /path/to/regent.py algorithm.rg

## Program Arguments ##
To control the thread counts with regent, the following arguments are fed to the program:
`-ll:cpu X` runs X compute threads
`-ll:util Y` runs Y utility threads
For best performance, we have found the X+Y equalling the number of cores,
and Y usually being 4 works well for single node runs.
`ll:csize Z` requests Z MB of CPU DRAM (per node) for the application, 
without this being set to a sufficient value the program will fail to run, 
though the default is sufficient for the 258x258 case.
There are also some specific parameters to tune the application:
`t1 A` and `t2 B` controls how the system is partitioned into work across the node.
We would recommend starting with `A=256` and `B = problemsize/2`, i.e.
    
    > /path/to/regent.py algorithm.rg -ll:cpu 8 -ll:cpu 4 -t1 256 -t2 512

## Output ##
The application output (if requested) is saved in `go2d_<time-step>_00000.dat`, 
and is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.

