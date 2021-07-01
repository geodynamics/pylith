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
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

subpackage = testing
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

TESTS = test_utils

check_PROGRAMS = test_utils

# Primary source files
test_utils_SOURCES = \
	TestStubMethodTracker.cc \
	test_driver.cc

noinst_HEADERS = \
	TestStubMethodTracker.hh

AM_CPPFLAGS += $(PETSC_CC_INCLUDES)
AM_CPPFLAGS += $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR)

test_utils_LDFLAGS = $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)
test_utils_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  test_utils_LDADD += -lnetcdf
endif


# End of file
