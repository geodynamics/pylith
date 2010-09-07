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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = topology
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	Distributor.hh \
	FieldBase.hh \
	Field.hh \
	Field.icc \
	Field.cc \
	Fields.hh \
	Fields.icc \
	FieldsNew.hh \
	FieldsNew.icc \
	FieldsNew.cc \
	Jacobian.hh \
	Mesh.hh \
	Mesh.icc \
	MeshOps.hh \
	RefineUniform.hh \
	ReverseCuthillMcKee.hh \
	SolutionFields.hh \
	SubMesh.hh \
	SubMesh.icc \
	SubMesh.cc \
	UniformSectionDS.hh \
	UniformSectionDS.cc \
	topologyfwd.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
