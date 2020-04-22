#!/usr/bin/env bash

# This script executes the selected matrix vector multiplication implementations
# with and without auto-vectorization enabled in the compiler.

output=vectorization.txt

echo "Vectorization plot" | tee $output

echo "Version VectorTime NoVecTime" >> $output
for v in `make versions`; do
    vec=$(OMP_NUM_THREADS=1 ./kdriver.avx2.$v | grep "Loop time" | awk '{print $5}')
    novec=$(OMP_NUM_THREADS=1 ./kdriver.novec.$v | grep "Loop time" | awk '{print $5}')
    echo $v $vec $novec >> $output
done
