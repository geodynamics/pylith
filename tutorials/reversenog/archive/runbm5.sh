#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith in tutorial using SCEC benchmark 5.
#
# ======================================================================

if [ $# != 1 ]; then
  echo "usage: runbm5.sh NPROCS"
  exit 1
fi
nprocs=$1

simroot="bm5"
dupext="fuldat prop statevar time"
sinext="coord connect bc split"

echo "Setting up symbolic links with prefix ${simroot}_${nprocs}..."
for ext in $sinext; do
  ln -s $simroot.$ext ${simroot}_$nprocs.$ext
done
for ext in $dupext; do
  for (( i=0; i < $nprocs; i+=1 )); do
    ln -s $simroot.$ext ${simroot}_$nprocs.$i.$ext
  done
done

echo "Running PyLith..."
cmd="mpirun -np $nprocs `which pylith3dapp.py` \
  --typos=relaxed \
  --scanner.fileRoot=${simroot}_$nprocs \
  --scanner.asciiOutput=full \
  --scanner.ucdOutput=ascii \
  -log_summary -pc_type bjacobi -sub_pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-09"
echo $cmd
eval $cmd

# end of file
