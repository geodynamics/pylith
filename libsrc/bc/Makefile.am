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
	DirichletBC.hh \
	DirichletBC.icc \
	DirichletBoundary.hh \
	DirichletBoundary.icc \
	Neumann.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
