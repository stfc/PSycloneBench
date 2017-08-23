#!/bin/bash
# Environment vars set-up in order to build dl_timer with the Intel compiler
export F90=ifort
export F90FLAGS="-fast -warn all"
export LDFLAGS=
export CC=icc
export CFLAGS="-std=c99"
export MPIF90=mpiifort
export OMPFLAGS="-openmp"

