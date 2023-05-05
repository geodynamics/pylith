# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
include $(top_srcdir)/tests/check_cppunit.am

SUBDIRS = data

TESTS = test_topology

check_PROGRAMS = test_topology

# Primary source files
test_topology_SOURCES = \
	TestMesh.cc \
	TestMeshOps.cc \
	TestSubmesh.cc \
	TestSubmesh_Cases.cc \
	TestFieldBase.cc \
	TestFieldMesh.cc \
	TestFieldMesh_Cases.cc \
	TestFieldQuery.cc \
	TestFieldQuery_Cases.cc \
	TestRefineUniform.cc \
	TestRefineUniform_Cases.cc \
	TestReverseCuthillMcKee.cc \
	TestReverseCuthillMcKee_Cases.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_cppunit.cc

#	TestFieldSubmesh.cc



dist_noinst_HEADERS = \
	TestMesh.hh \
	TestSubmesh.hh \
	TestMeshOps.hh \
	TestFieldBase.hh \
	TestFieldMesh.hh \
	TestFieldSubmesh.hh \
	TestFieldQuery.hh \
	TestRefineUniform.hh \
	TestReverseCuthillMcKee.hh



noinst_TMP = \
	jacobian.mat \
	jacobian.mat.info \
	mesh.tex \
	mesh.txt \
	mesh.vtk \
	mesh.vtu \
	mesh_petsc.h5 \
	mesh_xdmf.h5


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am



# End of file
