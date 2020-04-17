#!/usr/bin/env bash

output=vectorization.txt

echo "Vectorization plot" | tee $output

make && make kdriver.novec

echo "Version VectorTime NoVecTime" >> $output
for v in `make versions`; do
    #vec=$(OMP_NUM_THREADS=1 ./kdriver.avx512.$v | grep "Loop time" | awk '{print $5}')
    vec=$(OMP_NUM_THREADS=1 ./kdriver.avx2.$v | grep "Loop time" | awk '{print $5}')
    novec=$(OMP_NUM_THREADS=1 ./kdriver.novec.$v | grep "Loop time" | awk '{print $5}')
    echo $v $vec $novec >> $output
done
