#!/bin/bash

# Run several full-scale tests. Execute this script from the tests
# directory.

tests_dir=`pwd`

# Clean output directory
clean_output() {
  cd ${tests_dir}/${dir}
  if [ -d output ]; then rm output/*.vtk output/*.h5 output/*.xmf; fi
}

# ----------------------------------------------------------------------
# 2d/frictionslide
dir="2d/frictionslide"
clean_output

echo "RUNNING $dir tension"
pylith tension.cfg --nodes=1 >& tension_np1.log
pylith tension.cfg --nodes=2 >& tension_np2.log
pylith tension.cfg --nodes=3 >& tension_np3.log

echo "RUNNING $dir ratestate_weak"
pylith ratestate.cfg ratestate_weak.cfg --nodes=1 >& weak_np1.log
pylith ratestate.cfg ratestate_weak.cfg --nodes=2 >& weak_np2.log
pylith ratestate.cfg ratestate_weak.cfg --nodes=3 >& weak_np3.log

echo "RUNNING $dir ratestate_stable"
pylith ratestate.cfg ratestate_stable.cfg --nodes=1 >& stable_np1.log
pylith ratestate.cfg ratestate_stable.cfg --nodes=2 >& stable_np2.log
pylith ratestate.cfg ratestate_stable.cfg --nodes=3 >& stable_np3.log

# ----------------------------------------------------------------------
# 3d/cyclicfriction
dir="3d/cyclicfriction"
clean_output

echo "RUNNING $dir ASM"
pylith --nodes=1 >& asm_np1.log
pylith --nodes=2 >& asm_np2.log
pylith --nodes=5 >& asm_np5.log

echo "RUNNING $dir fieldsplit"
pylith fieldsplit.cfg --nodes=1 >& fieldsplit_np1.log
pylith fieldsplit.cfg --nodes=2 >& fieldsplit_np2.log
pylith fieldsplit.cfg --nodes=5 >& fielssplit_np5.log

# ----------------------------------------------------------------------
# Return to examples dir
cd ${tests_dir}
  
exit 0
