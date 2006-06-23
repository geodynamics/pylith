#!/bin/bash
# ======================================================================
#
# Shell script to run LithoMop using SCEC benchmark 1 (medium resolution).
#
# ======================================================================

if [ $# != 0 ]; then
  echo "usage: runbm1a-lithomop.sh"
  exit 1
fi
suff="lm"

simroot="bm1a"
dupext="fuldat prop statevar time"
sinext="coord connect bc split"

pyreflags="-typos=relaxed"
echo "Pyre flags:"
echo $pyreflags

lithomopflags="-lm3dscan.fileRoot=${simroot}_$suff \
  -lm3dscan.asciiOutput=none \
  -lm3dscan.ucdOutput=binary"
echo "LithoMop flags:"
echo $lithomopflags

petscflags="-log_summary \
  -pc_type bjacobi \
  -sub_pc_type ilu \
  -ksp_monitor \
  -ksp_rtol 1e-09"
echo "PETSc flags:"
echo $petscflags

echo "Setting up symbolic links with prefix ${simroot}_${suff}..."
for ext in $sinext; do
  ln -s $simroot.$ext ${simroot}_$suff.$ext
done
for ext in $dupext; do
  ln -s $simroot.$ext ${simroot}_$suff.$ext
done

echo "Running LithoMop..."
cmd="lithomop3dapp.py $pyreflags $lithomopflags $petscflags"

echo $cmd
eval $cmd

# end of file
