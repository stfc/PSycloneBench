# README #

This directory contains the 'multi-kernel' version of the tra_adv benchmark.
It is so named because the code has been broken-up into ten subroutines,
each containing a single stencil operation. This is currently necessary in
order for the Dawn/Dusk toolchain to generate correct, optimised CUDA code.

## Building ##

The code uses the [dl_timer](https://bitbucket.org/apeg/dl_timer)
library for timing. Provided you followed the steps in the top-level
[README](../../../../README.md#obtaining-the-code) then the source code
for this library will have been obtained when you cloned the
repository. The library is built automatically as part of the
compilation process.

The compiler and compiler flags to use must be set through the following
environment variables:

* F90      - the command with which to invoke the Fortran compiler
* F90FLAGS - flags to pass to the compiler, e.g. `-g`

Scripts to do this for various compilers may be found in the `compiler_setup`
directory at the root of this repository.

For the various targets that involve OpenACC, obviously an OpenACC-capable
compiler is required (such as the one from NVIDIA) and suitable flags
must be specified, see e.g. `../../../../compiler_setup/nvidia.sh`.

A Makefile is provided that defines various targets as follows:

|Target            | Description                                       |
|------------------|---------------------------------------------------|
|`tra_adv_serial`  | Serial version                                    |
|`tra_adv_acc`     | PSyclone-generated OpenACC version                |
|`tra_adv_acc_prof`| PSyclone-generated OpenACC with profiling[^1]     |
|`tra_adv_sir`     | Serial, PSyclone-transformed SIR-compliant version|
|`tra_adv_sir_acc` | PSyclone generated OpenACC, SIR-compliant version |

[^1]: This target uses the PSyclone profiling API and thus requires that a
  suitable wrapper library has been compiled (such as the one in
  PSyclone/lib/profiling/nvidia). The location of this library must be
  supplied via the `PSYCLONE_NVIDIA_LIB_DIR` environment variable.

Each of these targets will produce a directory of the same name containing
the transformed code and a compiled executable.

## Running ##

The problem size and number of iterations/time-steps to perform are 
set using environment variables:

* JPI - x-dimension of the domain
* JPJ - y-dimension of the domain
* JPK - number of vertical levels
* IT - number of iterations to perform

e.g.:

```
    cd tra_adv_serial
    JPI=130 JPJ=130 JPK=31 IT=100 ./tra_adv.exe
```