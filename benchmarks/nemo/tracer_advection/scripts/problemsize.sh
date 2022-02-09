#!/bin/env bash

# Bash script to execute the tracer-advection benchmark with increasing
# domain sizes.

if [ "$#" -lt 1 ] || [ ! -x "$1" ]; then
    echo "Wrong arguments. Usage: ../../problemsize.sh ./executable"
    exit
fi

echo "Running problem size in $PWD with $@" 
echo "N time/step" 

# Number of iterations to perform
export IT=500
# Number of vertical levels
export JPK=75

base=2
#for power in $(seq 4 12); do
for power in $(seq 4 9); do

    size=$(echo "$base^$power" | bc)
    export JPI=${size}
    export JPJ=${size}

    time=$(taskset -c 2 $@  | awk '{if ($1 == "Time-stepping") {print $5} }')

    echo $size $time
done
