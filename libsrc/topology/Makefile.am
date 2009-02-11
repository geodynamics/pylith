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
	Field.cc \
	Fields.hh \
	Fields.icc \
	Mesh.hh \
	Mesh.icc \
	MeshOps.hh \
	MeshRefiner.hh \
	RefineUniform.hh \
	SubMesh.hh \
	SubMesh.icc \
	SubMesh.cc \
	RestrictVisitor.hh \
	RestrictVisitor.icc \
	topologyfwd.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
