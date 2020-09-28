

if [ "$#" -lt 1 ] || [ ! -x "$1" ]; then
    echo "Wrong arguments. Usage: ../../mpi_weak_scalability.sh ./executable"
    exit
fi

echo "Running problem size in $PWD with $@" 
echo "N time/step Gb/s time/step/problemsize" 

line_i=4
line_j=5 

nthreads=32
base=2
base_size=2048

for power in $(seq 0 6); do

    nprocs=$(echo "$base^$power" | bc)

    # Alternate the size that grows, eg: 512x512, 1024x512, 1024x1024, 2048x1024, ...
    if (( $power%2 )); then
        size_i=$(echo "${base_size}*($base^($power/2 + 1))" | bc)
        size_j=$(echo "${base_size}*($base^($power/2))" | bc)
    else
        size_i=$(echo "${base_size}*($base^($power/2))" | bc)
        size_j=$(echo "${base_size}*($base^($power/2))" | bc)
    fi

    total_size=$(echo "${size_i}*${size_j}" | bc)
    sed -i --follow-symlinks "${line_i}s/.*/jpiglo = ${size_i}/" namelist
    sed -i --follow-symlinks "${line_j}s/.*/jpjglo = ${size_j}/" namelist
    
    export OMP_NUM_THREADS=$nthreads
    time=$(mpirun -n ${nprocs} -ppn 1 $@  | awk '{if ($1 == "Time-stepping") {print $5} }')

    echo $nprocs $size_i $size_j $total_size $time
done
