#!/bin/bash
echo "Let's see if this works "
pkill imaster
sleep 0.5
../master/imaster -v &
sleep 1
mpirun -n 4 ../slave/islave &
sleep 1
../jsub/jsub -v -cfg synth_cube.cfg
echo "All should be running now!"
echo "Do tail -f invlog.00001"
