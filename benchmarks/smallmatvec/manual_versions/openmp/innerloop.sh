#!/usr/bin/env bash

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
