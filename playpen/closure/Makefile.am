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

bin_PROGRAMS = testlagrit testcubit

testlagrit_SOURCES = \
	TestClosure.cc \
	test_lagrit.cc

testcubit_SOURCES = \
	TestClosure.cc \
	test_cubit.cc

noinst_HEADERS = \
	TestClosure.hh


AM_CPPFLAGS = \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testlagrit_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testlagrit_LDADD = \
	-ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

testcubit_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testcubit_LDADD = \
	-ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS) \
	-lnetcdf_c++ -lnetcdf


# End of file 
