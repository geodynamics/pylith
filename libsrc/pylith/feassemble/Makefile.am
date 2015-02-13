# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = feassemble
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	CellGeometry.hh \
	CellGeometry.icc \
	Constraint.hh \
	Constraint.icc \
	ElasticityExplicit.hh \
	ElasticityExplicitTri3.hh \
	ElasticityExplicitTet4.hh \
	ElasticityExplicitLgDeform.hh \
	ElasticityImplicit.hh \
	ElasticityImplicitLgDeform.hh \
	Integrator.hh \
	Integrator.icc \
	IntegratorElasticity.hh \
	IntegratorElasticityLgDeform.hh \
	GeometryLine2D.hh \
	GeometryLine3D.hh \
	GeometryTri2D.hh \
	GeometryTri3D.hh \
	GeometryTet3D.hh \
	GeometryQuad2D.hh \
	GeometryQuad3D.hh \
	GeometryHex3D.hh \
	QuadratureRefCell.hh \
	QuadratureRefCell.icc \
	QuadratureEngine.hh \
	QuadratureEngine.icc \
	Quadrature.hh \
	Quadrature.icc \
	Quadrature1Din2D.hh \
	Quadrature1Din2D.icc \
	Quadrature1Din3D.hh \
	Quadrature1Din3D.icc \
	Quadrature2D.hh \
	Quadrature2D.icc \
	Quadrature2Din3D.hh \
	Quadrature2Din3D.icc \
	Quadrature3D.hh \
	Quadrature3D.icc \
	feassemblefwd.hh

noinst_HEADERS =

dist_noinst_DATA = \
	jacobian2d_lgdeform.wxm \
	jacobian2d_nonsymm_lgdeform.wxm \
	jacobian3d_lgdeform.wxm \
	jacobian3d_nonsymm_lgdeform.wxm \
	tri3_elasticity.wxm \
	tet4_elasticity.wxm


# End of file 
