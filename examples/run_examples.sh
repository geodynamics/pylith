#!/bin/bash

# Run all of the examples. Execute this script from the examples
# directory.

examples_dir=`pwd`

# Run specified list of examples (single .cfg file)
run_examples() {
  nprocs=1
  if [ $# == 1 ]; then
    nprocs=$1
  fi
  cd ${examples_dir}/$dir
  rm *.vtk
  if [ -d output ]; then rm output/*.vtk; fi
  for example in $examples; do
    echo "RUNNING $dir/$example"
    pylith --nodes=$nprocs $example
  done
}

# Run specified list of examples (additional common .cfg files)
run_examples2() {
  commoncfg=$1
  cd ${examples_dir}/$dir
  for example in $examples; do
    echo "RUNNING $dir/$example"
    pylith $commoncfg $example
  done
}

# ----------------------------------------------------------------------
# twotri3
dir="twocells/twotri3"
examples="axialdisp.cfg sheardisp.cfg dislocation.cfg"
run_examples

# twoquad4
dir="twocells/twoquad4"
examples="axialdisp.cfg sheardisp.cfg axialtract.cfg dislocation.cfg"
run_examples

# twotet4
dir="twocells/twotet4"
examples="axialdisp.cfg dislocation.cfg"
run_examples

# twohex8
dir="twocells/twohex8"
examples="axialdisp.cfg sheardisp.cfg dislocation.cfg"
run_examples

# twotet4-geoproj
dir="twocells/twotet4-geoproj"
examples="dislocation.cfg"
run_examples

# ----------------------------------------------------------------------
# 3d/tet4
dir="3d/tet4"
examples="step01.cfg step02.cfg step03.cfg step04.cfg"
run_examples
run_examples 2
run_examples 3
run_examples 4
run_examples 5

# 3d/hex8
dir="3d/hex8"
examples="step01.cfg step02.cfg step03.cfg step04.cfg step05.cfg step06.cfg step07.cfg step08.cfg step09.cfg step10.cfg step11.cfg step12.cfg step13.cfg step14.cfg step15.cfg step16.cfg step17.cfg step18.cfg step19.cfg step20.cfg"
run_examples

examples="step01.cfg step03.cfg step06.cfg step15.cfg step19.cfg step20.cfg"
run_examples 2

# ----------------------------------------------------------------------
# subduction
dir="2d/subduction"
examples="step01.cfg step02.cfg step03.cfg step04.cfg"
run_examples
run_examples 2
run_examples 4
run_examples 5

# ----------------------------------------------------------------------
# bar_shearwave/tri3
dir="bar_shearwave/tri3"
examples="pylithapp.cfg"
run_examples

# bar_shearwave/quad4
dir="bar_shearwave/quad4"
examples="prescribedrup.cfg"
run_examples
examples="spontaneousrup_staticfriction.cfg spontaneousrup_slipweakening.cfg spontaneousrup_ratestateageing.cfg"
run_examples2 "spontaneousrup.cfg"

# bar_shearwave/tet4
dir="bar_shearwave/tet4"
examples="pylithapp.cfg"
run_examples

# bar_shearwave/hex8
dir="bar_shearwave/hex8"
examples="pylithapp.cfg"
run_examples


# ----------------------------------------------------------------------
# greensfns
dir="2d/greensfns/strikeslip"
examples="eqsim.cfg --problem=pylith.problems.GreensFns"
run_examples

dir="2d/greensfns/reverse"
examples="eqsim.cfg --problem=pylith.problems.GreensFns"
run_examples

# ----------------------------------------------------------------------
# Return to examples dir
cd ${examples_dir}
  
exit 0
