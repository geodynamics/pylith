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
	BoundaryCondition.hh \
	BoundaryCondition.icc \
	BoundaryConditionPoints.hh \
	BCIntegratorSubMesh.hh \
	TimeDependent.hh \
	TimeDependent.icc \
	TimeDependentPoints.hh \
	TimeDependentPoints.icc \
	AbsorbingDampers.hh \
	DirichletBC.hh \
	DirichletBC.icc \
	DirichletBoundary.hh \
	DirichletBoundary.icc \
	Neumann.hh \
	PointForce.hh \
	PointForce.icc \
	bcfwd.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
