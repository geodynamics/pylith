# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


include $(top_srcdir)/tests/check_catch2.am

SUBDIRS = data

TESTS = libtest_faults

check_PROGRAMS = libtest_faults

libtest_faults_SOURCES = \
	TestTransformTopology.cc \
	TestTransformTopology_Tri.cc \
	TestTransformTopology_Quad.cc \
	TestTransformTopology_Tet.cc \
	TestTransformTopology_Hex.cc \
	TestFaultCohesive.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


dist_noinst_HEADERS = \
	TestTransformTopology.hh


# End of file
