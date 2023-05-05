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

TESTS = test_bc

check_PROGRAMS = test_bc

# Primary source files
test_bc_SOURCES = \
	$(top_srcdir)/tests/src/driver_cppunit.cc

#	TestBoundaryMesh.cc \
#	TestBoundaryMesh_Cases.cc \
#	TestDirichletTimeDependent.cc \
#	TestDirichletTimeDependent_Cases.cc \
#	TestNeumannTimeDependent.cc \
#	TestNeumannTimeDependent_Cases.cc \
#	TestAbsorbingDampers.cc \
#	TestAbsorbingDampers_Cases.cc


dist_noinst_HEADERS = \
	TestBoundaryMesh.hh \
	TestDirichletTimeDependent.hh \
	TestNeumannTimeDependent.hh \
	TestAbsorbingDampers.hh


# End of file
