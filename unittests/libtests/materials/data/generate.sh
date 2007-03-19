#!/bin/bash
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

if (( $# != 1 )); then
  echo "usage: generate.sh elastic|viscoelastic|all"
  exit 1
fi

# //////////////////////////////////////////////////////////////////////
if [ $1 == "elastic" ] || [ $1 == "all" ]; then

  # 3-D ----------------------------------------------------------------

  python ElasticIsotropic3D.py \
    --data.namespace=pylith,materials \
    --data.object=ElasticIsotropic3DData \
    --data.parent=ElasticMaterialData

fi

# //////////////////////////////////////////////////////////////////////
if [ $1 == "viscoelastic" ] || [ $1 == "all" ]; then
  echo "" >& /dev/null

fi


# End of file 
