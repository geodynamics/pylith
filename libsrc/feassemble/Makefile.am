# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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
	ElasticityImplicit.hh \
	Integrator.hh \
	Integrator.icc \
	Elasticity.hh \
	GeometryPoint1D.hh \
	GeometryPoint2D.hh \
	GeometryPoint3D.hh \
	GeometryLine1D.hh \
	GeometryLine2D.hh \
	GeometryLine3D.hh \
	GeometryTri2D.hh \
	GeometryTri3D.hh \
	GeometryTet3D.hh \
	GeometryQuad2D.hh \
	GeometryQuad3D.hh \
	GeometryHex3D.hh \
	Quadrature.hh \
	Quadrature.icc \
	Quadrature0D.hh \
	Quadrature0D.icc \
	Quadrature1D.hh \
	Quadrature1D.icc \
	Quadrature1Din2D.hh \
	Quadrature1Din2D.icc \
	Quadrature1Din3D.hh \
	Quadrature1Din3D.icc \
	Quadrature2D.hh \
	Quadrature2D.icc \
	Quadrature2Din3D.hh \
	Quadrature2Din3D.icc \
	Quadrature3D.hh \
	Quadrature3D.icc

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
