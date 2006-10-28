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

  python Quadrature1DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1DLinear \
    --data.parent=QuadratureData

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
if [ $1 == "integrator" ] || [ $1 == "all" ]; then

  # 1-D ----------------------------------------------------------------

  python IntegratorInertia1DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=IntegratorDataInertia1DLinear \
    --data.parent=IntegratorData

  python IntegratorInertia1DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=IntegratorDataInertia1DQuadratic \
    --data.parent=IntegratorData

  # 2-D ----------------------------------------------------------------

  python IntegratorInertia2Din3DOne.py \
    --data.namespace=pylith,feassemble \
    --data.object=IntegratorDataInertia2Din3DOne \
    --data.parent=IntegratorData

  python IntegratorInertia2Din3DTwo.py \
    --data.namespace=pylith,feassemble \
    --data.object=IntegratorDataInertia2Din3DTwo \
    --data.parent=IntegratorData

  python IntegratorInertia2Din3DThree.py \
    --data.namespace=pylith,feassemble \
    --data.object=IntegratorDataInertia2Din3DThree \
    --data.parent=IntegratorData

  # 3-D ----------------------------------------------------------------

  python IntegratorInertia3DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=IntegratorDataInertia3DLinear \
    --data.parent=IntegratorData

  python IntegratorInertia3DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=IntegratorDataInertia3DQuadratic \
    --data.parent=IntegratorData

fi


# End of file 
