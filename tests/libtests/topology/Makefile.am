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
include $(top_srcdir)/tests/check_catch2.am

SUBDIRS = data

TESTS = libtest_topology

check_PROGRAMS = libtest_topology

# Primary source files
libtest_topology_SOURCES = \
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
	$(top_srcdir)/tests/src/driver_catch2.cc

#	TestFieldSubmesh.cc
#	TestFieldSubmesh.hh



dist_noinst_HEADERS = \
	TestSubmesh.hh \
	TestFieldMesh.hh \
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
