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
	GeometryPoint.hh \
	GeometryLine.hh \
	GeometryTri.hh \
	GeometryTet.hh \
	GeometryQuad.hh \
	GeometryHex.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
