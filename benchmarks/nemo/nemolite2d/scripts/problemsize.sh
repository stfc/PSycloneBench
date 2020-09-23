

if [ "$#" -lt 1 ] || [ ! -x "$1" ]; then
    echo "Wrong arguments. Usage: ../../problemsize.sh ./executable"
    exit
fi

echo "Running problem size in $PWD with $@" 
echo "N time/step Gb/s time/step/problemsize" 

line_i=4
line_j=5 

base=2
for power in $(seq 4 12); do

    size=$(echo "$base^$power" | bc)
    sed -i --follow-symlinks "${line_i}s/.*/jpiglo = ${size}/" namelist
    sed -i --follow-symlinks "${line_j}s/.*/jpjglo = ${size}/" namelist

    time=$($@  | awk '{if ($1 == "Time-stepping") {print $5} }')

    echo $size $time
done
