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

ACLOCAL_AMFLAGS = -I ./m4

bin_PROGRAMS = runtet4

# Primary source files
runtet4_SOURCES = \
	runtet4.cc

AM_CPPFLAGS = \
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

runtet4_LDADD = \
	-lpylith \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

CLEANFILES = $(bin_PROGRAMS)

# End of file 
