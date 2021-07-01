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

subpackage = bc
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = test_bc

check_PROGRAMS = test_bc

# Primary source files
test_bc_SOURCES = \
	test_driver.cc

#	TestBoundaryMesh.cc \
#	TestBoundaryMesh_Cases.cc \
#	TestDirichletTimeDependent.cc \
#	TestDirichletTimeDependent_Cases.cc \
#	TestNeumannTimeDependent.cc \
#	TestNeumannTimeDependent_Cases.cc \
#	TestAbsorbingDampers.cc \
#	TestAbsorbingDampers_Cases.cc


noinst_HEADERS = \
	TestBoundaryMesh.hh \
	TestDirichletTimeDependent.hh \
	TestNeumannTimeDependent.hh \
	TestAbsorbingDampers.hh


AM_CPPFLAGS += \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES)

test_bc_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

test_bc_LDADD = \
	-lcppunit \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  test_bc_LDADD += -lnetcdf
endif


# End of file
