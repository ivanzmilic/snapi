#!/bin/bash
echo "Running a simple synthesis according to the synth.cfg file..."
pkill imaster
sleep 0.5
../master/imaster -v &
echo "Master process running... "
sleep 0.5
../slave/islave &
echo "Worker process running..."
sleep 0.5
../jsub/jsub -v -cfg synth_3D1px.cfg
echo "Job submitted! Everything should be running now."
echo "Check the progress with tail -f invlog.000001."
echo "Once the calculation is finished, you can ../jsub/jub another job, or pkill imaster."