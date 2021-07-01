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

subpackage = topology
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = test_topology

check_PROGRAMS = test_topology

# Primary source files
test_topology_SOURCES = \
	TestMesh.cc \
	TestMeshOps.cc \
	TestSubmesh.cc \
	TestSubmesh_Cases.cc \
	TestFieldBase.cc \
	TestFieldMesh.cc \
	TestFieldMesh_Cases.cc \
	TestFieldQuery.cc \
	TestFieldQuery_Cases.cc \
	TestRefineUniform.cc \
	TestRefineUniform_Cases.cc \
	TestReverseCuthillMcKee.cc \
	TestReverseCuthillMcKee_Cases.cc \
	test_driver.cc

#	TestFieldSubmesh.cc



noinst_HEADERS = \
	TestMesh.hh \
	TestSubmesh.hh \
	TestMeshOps.hh \
	TestFieldBase.hh \
	TestFieldMesh.hh \
	TestFieldSubmesh.hh \
	TestFieldQuery.hh \
	TestRefineUniform.hh \
	TestReverseCuthillMcKee.hh


AM_CPPFLAGS += \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

test_topology_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  test_topology_LDADD += -lnetcdf
endif

noinst_tmp = \
	jacobian.mat \
	jacobian.mat.info

#CLEANFILES = $(noinst_tmp)


# End of file
