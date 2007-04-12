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

  # 2-D ----------------------------------------------------------------

  python ElasticPlaneStrain.py \
    --data.namespace=pylith,materials \
    --data.object=ElasticPlaneStrainData \
    --data.parent=ElasticMaterialData

  python ElasticPlaneStress.py \
    --data.namespace=pylith,materials \
    --data.object=ElasticPlaneStressData \
    --data.parent=ElasticMaterialData

  # 1-D ----------------------------------------------------------------

  python ElasticStrain1D.py \
    --data.namespace=pylith,materials \
    --data.object=ElasticStrain1DData \
    --data.parent=ElasticMaterialData

  python ElasticStress1D.py \
    --data.namespace=pylith,materials \
    --data.object=ElasticStress1DData \
    --data.parent=ElasticMaterialData

fi

# //////////////////////////////////////////////////////////////////////
if [ $1 == "viscoelastic" ] || [ $1 == "all" ]; then
  echo "" >& /dev/null

fi


# End of file 
