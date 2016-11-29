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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = utils
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

TESTS = testutils

check_PROGRAMS = testutils

# Primary source files
testutils_SOURCES = \
	TestEventLogger.cc \
	TestJournalingComponent.cc \
	TestPylithVersion.cc \
	TestPetscVersion.cc \
	TestDependenciesVersion.cc \
	test_utils.cc

noinst_HEADERS = \
	TestEventLogger.hh \
	TestJournalingComponent.hh \
	TestPylithVersion.hh \
	TestPetscVersion.hh \
	TestDependenciesVersion.hh

AM_CPPFLAGS += $(PETSC_CC_INCLUDES)
AM_CPPFLAGS += $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR)

testutils_LDFLAGS = $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)
testutils_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testutils_LDADD += -lnetcdf
endif


leakcheck: testutils
	valgrind --log-file=valgrind_utils.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/testutils


# End of file
