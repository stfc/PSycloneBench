#!/usr/bin/env bash

# This scripts executes serially (with 1 thread) the selected implementation
# with increasing problem sizes. The first loop changes the horizontal
# dimension and writes the results at 'problemsize.txt` the second loop
# increases the vertical dimension and writes the results at 'nlayers.txt'

# Update the for statements to change the target implementation and problem
# size ranges.

output=problemsize.txt

echo "Problem size plot" | tee $output

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
