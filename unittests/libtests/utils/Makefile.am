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
# Copyright (c) 2010-2015 University of California, Davis
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
	test_utils.cc

noinst_HEADERS = \
	TestEventLogger.hh

AM_CPPFLAGS += $(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES)

testutils_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testutils_LDADD += -lnetcdf_c++ -lnetcdf
endif


leakcheck: testutils
	valgrind --log-file=valgrind_utils.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/testutils


# End of file 
