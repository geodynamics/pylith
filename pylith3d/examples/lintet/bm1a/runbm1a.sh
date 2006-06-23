#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith in tutorial using SCEC benchmark 5.
#
# ======================================================================

if [ $# != 1 ]; then
  echo "usage: runbm1a.sh NPROCS"
  exit 1
fi
nprocs=$1

simroot="bm1a"
dupext="fuldat prop statevar time"
sinext="coord connect bc split"

pyreflags="-typos=relaxed"
echo "Pyre flags:"
echo $pyreflags

pylithflags="-pl3dscan.fileRoot=${simroot}_$nprocs \
  -pl3dscan.asciiOutput=none \
  -pl3dscan.ucdOutput=binary"
# Do not use pythonTimestep for now until all the bugs are worked out.
#  -pl3dscan.ucdOutput=ascii \
#  -pl3dscan.pythonTimestep=1 "
echo "PyLith flags:"
echo $pylithflags

petscflags="-log_summary \
  -pc_type bjacobi \
  -sub_pc_type ilu \
  -ksp_monitor \
  -ksp_rtol 1e-09"
echo "PETSc flags:"
echo $petscflags

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
cmd="mpiexec -np $nprocs pylith3dapp.py $pyreflags $pylithflags $petscflags"

echo $cmd
eval $cmd

# end of file
