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

subpackage = bc
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	AbsorbingDampers.hh \
	BoundaryCondition.hh \
	BoundaryCondition.icc \
	DirichletBoundary.hh \
	DirichletBoundary.icc \
	DirichletPoints.hh \
	DirichletPoints.icc \
	Neumann.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
