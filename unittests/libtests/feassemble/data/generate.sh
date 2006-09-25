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

# 1-D ------------------------------------------------------------------

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

# 2-D ------------------------------------------------------------------

python Quadrature2DLinear.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2DLinear \
    --data.parent=QuadratureData

python Quadrature2DQuadratic.py \
    --data.namespace=pylith,feassemble \
    --data.object=QuadratureData2DQuadratic \
    --data.parent=QuadratureData

# End of file 
