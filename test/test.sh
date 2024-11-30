#!/bin/bash
echo "Let's see if this works "
pkill imaster
../master/imaster -v &
../slave/islave &
../jsub/jsub -v -cfg synth.cfg
echo "Done!"
