Build Configuration Files
=========================

The directory containing this file contains various bash scripts
that set-up the necessary environment variables used by the
Makefiles for the benchmarks. Typically these are:

F90 - the Fortran compiler to use
F90FLAGS - the flags to pass to the Fortran compiler
OMPFLAGS - the compiler flags to use to enable OpenMP
LDFLAGS - flags to pass to the linker (if any)
AR - the command to use when creating an archive (.a) file
     from object (.o) files
ARFLAGS - the flags to pass to the archiver.

The latter two are really just to support the Intel compiler since
we must use Intel's xiar when doing Inter-Procedural Optimisation.