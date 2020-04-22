#!/usr/bin/env bash

# This scripts simulates the performance of fusing the inner loop of
# multiple matrix multiplications together by repeating the multiplication
# loop on the same data an arbitrary number of times (defined by INNERREPS).
# Note that this only works with the 'nlayersf_moreops' version.

output=innerloop.txt

echo "innerloop plot" | tee $output

echo "Version VectorTime NoVecTime" >> $output
for n in 1 2 4 8 12 16 20 24 28 32 64 128; do
    ifort -qopenmp -Ofast -xHost -no-vec -no-simd -DINNERREPS=$n -qopt_report=5 -c matrix_vector_kernel_mod.F90 -o matrix_vector_kernel_mod.novec.o
    ifort -qopenmp -Ofast -xHost -DINNERREPS=$n -qopt_report=5 -c matrix_vector_kernel_mod.F90 -o matrix_vector_kernel_mod.o
    make && make kdriver.novec
    vec=$(OMP_NUM_THREADS=1 ./kdriver.avx512.nlayersf_moreops | grep "Loop time" | awk '{print $5}')
    novec=$(OMP_NUM_THREADS=1 ./kdriver.novec.nlayersf_moreops | grep "Loop time" | awk '{print $5}')
    echo $n $vec $novec >> $output
done
