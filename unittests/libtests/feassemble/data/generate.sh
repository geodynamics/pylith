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

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1DQuadratic \
    --data.parent=QuadratureData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din2DLinear \
    --data.parent=QuadratureData \
    --mesh=Mesh1Din2DLinear \
    --quadrature=Quadrature1Din2DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din2DQuadratic \
    --data.parent=QuadratureData \
    --mesh=Mesh1Din2DQuadratic \
    --quadrature=Quadrature1Din2DQuadratic

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din3DLinear \
    --data.parent=QuadratureData \
    --mesh=Mesh1Din3DLinear \
    --quadrature=Quadrature1Din3DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData1Din3DQuadratic \
    --data.parent=QuadratureData \
    --mesh=Mesh1Din3DQuadratic \
    --quadrature=Quadrature1Din3DQuadratic

  # 2-D ----------------------------------------------------------------

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2DLinear \
    --data.parent=QuadratureData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2DQuadratic \
    --data.parent=QuadratureData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearXYZ \
    --data.parent=QuadratureData \
    --mesh=Mesh2Din3DLinearXYZ \
    --quadrature=Quadrature2Din3DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearXY \
    --data.parent=QuadratureData \
    --mesh=Mesh2Din3DLinearXY \
    --quadrature=Quadrature2Din3DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearYZ \
    --data.parent=QuadratureData \
    --mesh=Mesh2Din3DLinearYZ \
    --quadrature=Quadrature2Din3DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DLinearXZ \
    --data.parent=QuadratureData \
    --mesh=Mesh2Din3DLinearXZ \
    --quadrature=Quadrature2Din3DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2Din3DQuadratic \
    --data.parent=QuadratureData \
    --mesh=Mesh2Din3DQuadratic \
    --quadrature=Quadrature2Din3DQuadratic

  # 3-D ----------------------------------------------------------------

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData3DLinear \
    --data.parent=QuadratureData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear

  python QuadratureApp.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData3DQuadratic \
    --data.parent=QuadratureData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic

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
