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

TESTS = libtest_feassemble

check_PROGRAMS = libtest_feassemble

libtest_feassemble_SOURCES = \
	TestAuxiliaryFactory.cc \
	TestInterfacePatches.cc \
	TestInterfacePatches_Quad.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


dist_noinst_HEADERS = \
	TestInterfacePatches.hh


# End of file
