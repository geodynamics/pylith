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

TESTS = libtest_testing

check_PROGRAMS = libtest_testing

# Primary source files
libtest_testing_SOURCES = \
	TestStubMethodTracker.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc

dist_noinst_HEADERS = 

# End of file
