

if [ "$#" -lt 1 ] || [ ! -x "$1" ]; then
    echo "Wrong arguments. Usage: ../../mpi_weak_scalability.sh ./executable"
    exit
fi

echo "Running problem size in $PWD with $@" 
echo "N time/step Gb/s time/step/problemsize" 

line_i=4
line_j=5 

nthreads=2
ranks_per_node=$(echo "32/$nthreads" | bc)
base=2
base_size=2048

for power in $(seq 0 8); do

    nprocs=$(echo "($base^$power)*$ranks_per_node" | bc)

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
    export I_MPI_PIN_DOMAIN=omp
    time=$(mpirun -n ${nprocs} -ppn ${ranks_per_node} $@  | awk '{if ($1 == "Time-stepping") {print $5} }')

    echo $nprocs $size_i $size_j $total_size $time
done
