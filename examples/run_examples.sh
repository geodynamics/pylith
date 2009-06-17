#!/bin/bash

# Run all of the examples. Execute this script from the examples
# directory.

examples_dir=`pwd`

# Run specified list of examples
run_examples() {
  cd ${examples_dir}/$dir
  rm *.vtk
  if [ -d output ]; then rm output/*.vtk; fi
  for example in $examples; do
    echo "RUNNING $dir/$example"
    pylith $example
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
examples="shearxy.cfg dislocation.cfg"
run_examples

# 3d/hex8
dir="3d/hex8"
examples="shearxy.cfg dislocation.cfg gravity.cfg gravity_istress.cfg savageprescott.cfg"
run_examples

# ----------------------------------------------------------------------
# bar_shearwave/tri3
dir="bar_shearwave/tri3"
examples="pylithapp.cfg"
run_examples

# bar_shearwave/quad4
dir="bar_shearwave/quad4"
examples="pylithapp.cfg"
run_examples

# bar_shearwave/tet4
dir="bar_shearwave/tet4"
examples="pylithapp.cfg"
run_examples

# bar_shearwave/hex8
dir="bar_shearwave/hex8"
examples="pylithapp.cfg"
run_examples



# Return to examples dir
cd ${examples_dir}
  
exit 0
