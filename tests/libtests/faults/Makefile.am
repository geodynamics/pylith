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

TESTS = libtest_faults

check_PROGRAMS = libtest_faults

libtest_faults_SOURCES = \
	TestAdjustTopology.cc \
	TestAdjustTopology_Tri.cc \
	TestAdjustTopology_Quad.cc \
	TestAdjustTopology_Tet.cc \
	TestAdjustTopology_Hex.cc \
	TestFaultCohesive.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


dist_noinst_HEADERS = \
	TestAdjustTopology.hh


# End of file
