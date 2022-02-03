# README #

This directory contains the 'compute-in-subroutine' version of the
tra_adv benchmark. It is so named because the original, single-file
version of the benchmark has been broken-up into two routines - the
top-level driver program and a subroutine containing all computation.
This is required in order to construct the Stencil Intermediate
Representation (used by the Dawn/Dusk toolchain) of the code since
that IR only supports compute constructs.

Note that limitations in the SIR produced by PSyclone and in the Dawn
frontend mean that currently, correct CUDA code cannot be generated
for this version of the benchmark.

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

|Target          | Description                                                 |
|----------------|-------------------------------------------------------------|
|`tra_adv_serial`  | Serial Fortran version                                    |
|`tra_adv_no_auto_serial` | Serial with automatic work arrays moved to module scope |
|`tra_adv_acc`     | PSyclone-generated OpenACC Fortran version                |
|`tra_adv_acc_prof`| PSyclone-generated OpenACC Fortran with profiling[^1]     |
|`tra_adv_sir`     | Serial, PSyclone-transformed SIR-compliant Fortran version|
|`tra_adv_sir_acc` | PSyclone generated OpenACC, SIR-compliant Fortran version |
|`tra_adv_sir_acc_prof`| PSyclone generated OpenACC, SIR-compliant Fortran version with profiling[^1] |

[^1]: This target uses the PSyclone profiling API and thus requires that a
  suitable wrapper library has been compiled (such as the one in
  PSyclone/lib/profiling/nvidia). The location of this library must be
  supplied via the `PSYCLONE_NVIDIA_LIB_DIR` environment variable.

Each of these targets will produce a directory of the same name containing
the transformed code and a compiled executable.

As implemented, the 'compute' subroutine in `tra_adv_compute_auto_arrays.F90`
uses Fortran `automatic` arrays. It is known
(https://github.com/stfc/PSycloneBench/pull/76#issuecomment-943485412) that
this causes performance problems when using ManagedMemory on GPU. A potential
solution to this is shown in `tra_adv_compute.F90` where these arrays have
been promoted to module scope. This is the code that is built by the
`tra_adv_no_auto_serial` target.

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
