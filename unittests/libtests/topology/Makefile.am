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

subpackage = topology
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testtopology

check_PROGRAMS = testtopology

# Primary source files
testtopology_SOURCES = \
	TestMesh.cc \
	TestMeshOps.cc \
	TestSubMesh.cc \
	TestFieldBase.cc \
	TestFieldMesh.cc \
	TestFieldSubMesh.cc \
	TestFieldsMesh.cc \
	TestFieldsSubMesh.cc \
	TestSolutionFields.cc \
	TestJacobian.cc \
	TestRefineUniform.cc \
	TestReverseCuthillMcKee.cc \
	test_topology.cc



noinst_HEADERS = \
	TestMesh.hh \
	TestSubMesh.hh \
	TestMeshOps.hh \
	TestFieldBase.hh \
	TestFieldMesh.hh \
	TestFieldSubMesh.hh \
	TestFieldsMesh.hh \
	TestFieldsSubMesh.hh \
	TestSolutionFields.hh \
	TestRefineUniform.hh \
	TestReverseCuthillMcKee.hh \
	TestJacobian.hh



# Source files associated with testing data
testtopology_SOURCES += \
	data/MeshDataCohesive.cc \
	data/MeshDataCohesiveTri3Level1.cc \
	data/MeshDataCohesiveTri3Level1Fault1.cc \
	data/MeshDataCohesiveQuad4Level1.cc \
	data/MeshDataCohesiveQuad4Level1Fault1.cc \
	data/MeshDataCohesiveTet4Level1.cc \
	data/MeshDataCohesiveTet4Level1Fault1.cc \
	data/MeshDataCohesiveHex8Level1.cc \
	data/MeshDataCohesiveHex8Level1Fault1.cc

noinst_HEADERS += \
	data/MeshDataCohesive.hh \
	data/MeshDataCohesiveTri3Level1.hh \
	data/MeshDataCohesiveTri3Level1Fault1.hh \
	data/MeshDataCohesiveQuad4Level1.hh \
	data/MeshDataCohesiveQuad4Level1Fault1.hh \
	data/MeshDataCohesiveTet4Level1.hh \
	data/MeshDataCohesiveTet4Level1Fault1.hh \
	data/MeshDataCohesiveHex8Level1.hh \
	data/MeshDataCohesiveHex8Level1Fault1.hh

AM_CPPFLAGS += \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testtopology_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testtopology_LDADD += -lnetcdf_c++ -lnetcdf
endif

noinst_tmp = \
	jacobian.mat \
	jacobian.mat.info

CLEANFILES = $(noinst_tmp)


leakcheck: testtopology
	valgrind --log-file=valgrind_topology.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/testtopology


# End of file 
