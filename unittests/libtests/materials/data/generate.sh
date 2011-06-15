#!/bin/bash
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
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

  python MaxwellIsotropic3DElastic.py \
    --data.namespace=pylith,materials \
    --data.object=MaxwellIsotropic3DElasticData \
    --data.parent=ElasticMaterialData

  python PowerLaw3DElastic.py \
    --data.namespace=pylith,materials \
    --data.object=PowerLaw3DElasticData \
    --data.parent=ElasticMaterialData

  python GenMaxwellIsotropic3DElastic.py \
    --data.namespace=pylith,materials \
    --data.object=GenMaxwellIsotropic3DElasticData \
    --data.parent=ElasticMaterialData

  python GenMaxwellQpQsIsotropic3DElastic.py \
    --data.namespace=pylith,materials \
    --data.object=GenMaxwellQpQsIsotropic3DElasticData \
    --data.parent=ElasticMaterialData

  python DruckerPrager3DElastic.py \
    --data.namespace=pylith,materials \
    --data.object=DruckerPrager3DElasticData \
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

  python MaxwellPlaneStrainElastic.py \
    --data.namespace=pylith,materials \
    --data.object=MaxwellPlaneStrainElasticData \
    --data.parent=ElasticMaterialData

  python GenMaxwellPlaneStrainElastic.py \
    --data.namespace=pylith,materials \
    --data.object=GenMaxwellPlaneStrainElasticData \
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

  python GenMaxwellIsotropic3DTimeDep.py \
    --data.namespace=pylith,materials \
    --data.object=GenMaxwellIsotropic3DTimeDepData \
    --data.parent=ElasticMaterialData

  python GenMaxwellQpQsIsotropic3DTimeDep.py \
    --data.namespace=pylith,materials \
    --data.object=GenMaxwellQpQsIsotropic3DTimeDepData \
    --data.parent=ElasticMaterialData

  python MaxwellIsotropic3DTimeDep.py \
    --data.namespace=pylith,materials \
    --data.object=MaxwellIsotropic3DTimeDepData \
    --data.parent=ElasticMaterialData

  python GenMaxwellPlaneStrainTimeDep.py \
    --data.namespace=pylith,materials \
    --data.object=GenMaxwellPlaneStrainTimeDepData \
    --data.parent=ElasticMaterialData

  python MaxwellPlaneStrainTimeDep.py \
    --data.namespace=pylith,materials \
    --data.object=MaxwellPlaneStrainTimeDepData \
    --data.parent=ElasticMaterialData

  python PowerLaw3DTimeDep.py \
    --data.namespace=pylith,materials \
    --data.object=PowerLaw3DTimeDepData \
    --data.parent=ElasticMaterialData

  python DruckerPrager3DTimeDep.py \
    --data.namespace=pylith,materials \
    --data.object=DruckerPrager3DTimeDepData \
    --data.parent=ElasticMaterialData

fi


# End of file 
