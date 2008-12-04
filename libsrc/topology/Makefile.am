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
	Distributor.hh \
	Field.hh \
	Field.icc \
	FieldsManager.hh \
	FieldOps.hh \
	Mesh.hh \
	Mesh.icc \
	MeshOps.hh \
	MeshRefiner.hh \
	RefineUniform.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
