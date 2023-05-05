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

TESTS = test_utils

check_PROGRAMS = test_utils

# Primary source files
test_utils_SOURCES = \
	TestEventLogger.cc \
	TestPyreComponent.cc \
	TestGenericComponent.cc \
	TestPylithVersion.cc \
	TestPetscVersion.cc \
	TestDependenciesVersion.cc \
	$(top_srcdir)/tests/src/driver_cppunit.cc

dist_noinst_HEADERS = \
	TestEventLogger.hh \
	TestPyreComponent.hh \
	TestGenericComponent.hh \
	TestPylithVersion.hh \
	TestPetscVersion.hh \
	TestDependenciesVersion.hh


# End of file
