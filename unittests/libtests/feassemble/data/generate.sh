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
  echo "usage: generate.sh quadrature|integrator|all"
  exit 1
fi

# //////////////////////////////////////////////////////////////////////
if [ $1 == "quadrature" ] || [ $1 == "all" ]; then

  # 1-D ----------------------------------------------------------------

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1DLinear \
    --data.parent=QuadratureData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear

  python Quadrature1DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1DQuadratic \
    --data.parent=QuadratureData

  python Quadrature1Din2DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din2DLinear \
    --data.parent=QuadratureData

  python Quadrature1Din2DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din2DQuadratic \
    --data.parent=QuadratureData

  python Quadrature1Din3DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din3DLinear \
    --data.parent=QuadratureData

  python Quadrature1Din3DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din3DQuadratic \
    --data.parent=QuadratureData

  # 2-D ----------------------------------------------------------------

  python Quadrature2DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2DLinear \
    --data.parent=QuadratureData

  python Quadrature2DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2DQuadratic \
    --data.parent=QuadratureData

  python Quadrature2Din3DLinearXYZ.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearXYZ \
    --data.parent=QuadratureData

  python Quadrature2Din3DLinearXY.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearXY \
    --data.parent=QuadratureData

  python Quadrature2Din3DLinearYZ.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearYZ \
    --data.parent=QuadratureData

  python Quadrature2Din3DLinearXZ.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearXZ \
    --data.parent=QuadratureData

  python Quadrature2Din3DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DQuadratic \
    --data.parent=QuadratureData

  # 3-D ----------------------------------------------------------------

  python Quadrature3DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData3DLinear \
    --data.parent=QuadratureData

  python Quadrature3DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData3DQuadratic \
    --data.parent=QuadratureData

fi

# //////////////////////////////////////////////////////////////////////
if [ $1 == "explicit" ] || [ $1 == "all" ]; then

  # 1-D ----------------------------------------------------------------

  python ElasticityExplicit.py \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitData1DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  # 2-D ----------------------------------------------------------------


  # 3-D ----------------------------------------------------------------

fi


# End of file 
