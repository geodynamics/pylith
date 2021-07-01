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

subpackage = problems
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = test_problems

check_PROGRAMS = test_problems

# Primary source files
test_problems_SOURCES = \
	TestObserversSoln.cc \
	TestObserversPhysics.cc \
	TestSolutionFactory.cc \
	TestSolutionFactory_Cases.cc \
	TestProgressMonitor.cc \
	TestProgressMonitorTime.cc \
	test_driver.cc



noinst_HEADERS = \
	TestPhysics.hh \
	TestObserversSoln.hh \
	TestObserversPhysics.hh \
	TestSolutionFactory.hh \
	TestProgressMonitor.hh \
	TestProgressMonitor.hh

AM_CPPFLAGS += \
	$(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

test_problems_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

test_problems_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  test_problems_LDADD += -lnetcdf
endif

#CLEANFILES = $(noinst_tmp)


# End of file
