#!/usr/bin/env bash

output=problemsize.txt

echo "Problem size plot" | tee $output

#make && make kdriver.novec

#for v in `make versions`; do
for v in nlayersf ; do
    echo $v >> $output
    for size in 1 `seq 4 4 24`; do
        nlayers=40
        echo $v $size
        t=$(OMP_NUM_THREADS=1 ./kdriver.avx512.$v -gen $size $nlayers | grep "Loop time" | awk '{print $5}')
        echo $size $nlayers $t >> $output
    done
done

output=nlayers.txt

echo "Problem size plot" | tee $output

#for v in `make versions`; do
for v in nlayersf ; do
    echo $v >> $output
    for nlayers in 1 `seq 10 10 120`; do
        size=12
        echo $v $size
        t=$(OMP_NUM_THREADS=1 ./kdriver.avx512.$v -gen $size $nlayers | grep "Loop time" | awk '{print $5}')
        echo $size $nlayers $t >> $output
    done
done
