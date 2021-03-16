# Manual SYCL version of PSyKAl NEMOLite2D #

This directory contains a version of the PSyKAl form of NEMOLite2D
that has had its PSy layer manually parallelised using SYCL.
The Fortran(NemoLite2d) to C++/SYCL(PSy-layer) calls are done with the
interface introduced in the `psykal_cpp` version.

SYCL is a Khronos Standard for Performance Portable C++, several
implementations of this standard are available. This code aims to be
compliant with the standard and usable from any SYCL implementation but
it has been developed using the Intel OneAPI implementation and the correct
compilation with other implementations is not guaranteed.

## Compiling and Implementations ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. `../../../compiler_setup/gnu.sh` for sample
settings for the Gnu compiler suite. Additionally, it uses the `SYCL_COMP`
(defaults to `dpcpp`) and `SYCL_FLAGS` (defaults to `-lsycl -O3 -mtune=native`)
environment variables to select and tune the SYCL implementation to use.
The default values are for the Intel OneAPI implementation.

Note that the linking step is also done with the SYCL compiler, but as it links
also Fortran object files, the libgfortran library needs to be in the library path.

If you are using Bash and the GNU family of compilers then:

    > source ../../../../../compiler_setup/gnu.sh
    > make

will compile the default `nemolite2d_sycl.exe` binary.

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

SYCL can run on several target platforms (including accelerated devices),
by default it will run on the host platform. For example:

    > ./nemolite2d_sycl.exe

The execution will also list all the available SYCL platforms in the system,
to select a different target platform set the `SYCL_PLATFORM` environment
variable to the platform of choice, for example:

    > SYCL_PLATFORM=2 ./nemolite2d_sycl.exe

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.
