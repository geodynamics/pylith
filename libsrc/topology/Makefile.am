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

subpackage = topology
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	CellGeometry.hh \
	CellGeometry.icc \
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
	GeometryHex3D.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
