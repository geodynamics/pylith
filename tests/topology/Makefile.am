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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

if ENABLE_CUBIT

TESTS = testtopology

check_PROGRAMS = testtopology

# Primary source files
testtopology_SOURCES = \
	test_topology.cc


AM_CPPFLAGS = -I$(top_builddir)/include

AM_CPPFLAGS += \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testtopology_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testtopology_LDADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	-lnetcdf_c++ -lnetcdf \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)


endif


# End of file 
