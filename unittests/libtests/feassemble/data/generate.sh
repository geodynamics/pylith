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
  echo "usage: generate.sh quadrature|explicit|implicit|all"
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

  python ElasticityExplicitApp.py \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitData1DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityExplicitApp.py \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitData1DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic

    python ElasticityExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitGravData1DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitGravData1DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic

  # 2-D ----------------------------------------------------------------

  python ElasticityExplicitApp.py \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitData2DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityExplicitApp.py \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitData2DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic

  python ElasticityExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitGravData2DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitGravData2DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic



  # 3-D ----------------------------------------------------------------

  python ElasticityExplicitApp.py \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitData3DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityExplicitApp.py \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitData3DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic

  python ElasticityExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitGravData3DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitGravData3DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic
fi


# //////////////////////////////////////////////////////////////////////
if [ $1 == "implicit" ] || [ $1 == "all" ]; then

  # 1-D ----------------------------------------------------------------

  python ElasticityApp.py \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitData1DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityApp.py \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitData1DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic

  python ElasticityApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitGravData1DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitGravData1DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic

  # 2-D ----------------------------------------------------------------

  python ElasticityApp.py \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitData2DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityApp.py \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitData2DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic

  python ElasticityApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitGravData2DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitGravData2DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic


  # 3-D ----------------------------------------------------------------

  python ElasticityApp.py \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitData3DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityApp.py \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitData3DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic

  python ElasticityApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitGravData3DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicit \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitGravData3DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic

fi

# //////////////////////////////////////////////////////////////////////
if [ $1 == "explicitlgdeform" ] || [ $1 == "all" ]; then

  # 1-D ----------------------------------------------------------------

  python ElasticityLgDeformExplicitApp.py \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformData1DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityLgDeformExplicitApp.py \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformData1DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic

  python ElasticityLgDeformExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformGravData1DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityLgDeformExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformGravData1DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic


  # 2-D ----------------------------------------------------------------

  python ElasticityLgDeformExplicitApp.py \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformData2DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityLgDeformExplicitApp.py \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformData2DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic

  python ElasticityLgDeformExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformGravData2DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityLgDeformExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformGravData2DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic


  # 3-D ----------------------------------------------------------------

  python ElasticityLgDeformExplicitApp.py \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformData3DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityLgDeformExplicitApp.py \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformData3DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic

  python ElasticityLgDeformExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformGravData3DLinear \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityLgDeformExplicitApp.py \
    --use_gravity=True \
    --formulation=ElasticityExplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityExplicitLgDeformGravData3DQuadratic \
    --data.parent=ElasticityExplicitData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic

fi


# //////////////////////////////////////////////////////////////////////
if [ $1 == "implicitlgdeform" ] || [ $1 == "all" ]; then

  # 1-D ----------------------------------------------------------------

  python ElasticityLgDeformApp.py \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformData1DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityLgDeformApp.py \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformData1DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic

  python ElasticityLgDeformApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformGravData1DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh1DLinear \
    --quadrature=Quadrature1DLinear \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DLinear

  python ElasticityLgDeformApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformGravData1DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh1DQuadratic \
    --quadrature=Quadrature1DQuadratic \
    --material=MaterialElasticStrain1D \
    --solution=Solution1DQuadratic


  # 2-D ----------------------------------------------------------------

  python ElasticityLgDeformApp.py \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformData2DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityLgDeformApp.py \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformData2DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic

  python ElasticityLgDeformApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformGravData2DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh2DLinear \
    --quadrature=Quadrature2DLinear \
    --material=ElasticPlaneStrain \
    --solution=Solution2DLinear

  python ElasticityLgDeformApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformGravData2DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh2DQuadratic \
    --quadrature=Quadrature2DQuadratic \
    --material=ElasticPlaneStrain \
    --solution=Solution2DQuadratic


  # 3-D ----------------------------------------------------------------

  python ElasticityLgDeformApp.py \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformData3DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityLgDeformApp.py \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformData3DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic

  python ElasticityLgDeformApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformGravData3DLinear \
    --data.parent=IntegratorData \
    --mesh=Mesh3DLinear \
    --quadrature=Quadrature3DLinear \
    --material=ElasticIsotropic3D \
    --solution=Solution3DLinear

  python ElasticityLgDeformApp.py \
    --use_gravity=True \
    --formulation=ElasticityImplicitLgDeform \
    --data.namespace=pylith,feassemble \
    --data.object=ElasticityImplicitLgDeformGravData3DQuadratic \
    --data.parent=IntegratorData \
    --mesh=Mesh3DQuadratic \
    --quadrature=Quadrature3DQuadratic \
    --material=ElasticIsotropic3D \
    --solution=Solution3DQuadratic

fi


# End of file 
