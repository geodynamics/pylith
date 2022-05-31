#!/bin/bash

subdirs=("box-2d" "reverse-2d" "strikeslip-2d")

#for d in "box-2d box-3d reverse-2d strikeslip-2d subduction-2d subduction-3d": do
for d in ${subdirs[*]}; do
    pushd $d && pylith_cfgsearch --display=features --output-format=markdown --path=../../../../examples/$d && popd
done

# Remove superfluous metadata files (not used).
rm strikeslip-2d/step01_slip_cubit-synopsis.md