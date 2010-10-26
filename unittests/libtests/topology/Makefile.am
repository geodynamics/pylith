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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = topology
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testtopology

check_PROGRAMS = testtopology

# Primary source files
testtopology_SOURCES = \
	TestMesh.cc \
	TestMeshOps.cc \
	TestSubMesh.cc \
	TestFieldMesh.cc \
	TestFieldSubMesh.cc \
	TestFieldsMesh.cc \
	TestFieldsSubMesh.cc \
	TestSolutionFields.cc \
	TestJacobian.cc \
	test_topology.cc


noinst_HEADERS = \
	TestMesh.hh \
	TestSubMesh.hh \
	TestMeshOps.hh \
	TestFieldMesh.hh \
	TestFieldSubMesh.hh \
	TestFieldsMesh.hh \
	TestFieldsSubMesh.hh \
	TestSolutionFields.hh \
	TestJacobian.hh


# Source files associated with testing data
testtopology_SOURCES += 

noinst_HEADERS += 

AM_CPPFLAGS = \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testtopology_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testtopology_LDADD += -lnetcdf_c++ -lnetcdf
endif

noinst_tmp = \
	jacobian.mat \
	jacobian.mat.info

CLEANFILES = $(noinst_tmp)


# End of file 
