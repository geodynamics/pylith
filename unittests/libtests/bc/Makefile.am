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

include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testbc

check_PROGRAMS = testbc

# Primary source files
testbc_SOURCES = \
	TestBoundaryCondition.cc \
	TestBoundaryConditionPoints.cc \
	TestTimeDependent.cc \
	TestTimeDependentPoints.cc \
	TestBoundaryMesh.cc \
	TestBoundaryMeshCases.cc \
	TestAbsorbingDampers.cc \
	TestAbsorbingDampersCases.cc \
	TestDirichletBC.cc \
	TestDirichletBCCases.cc \
	TestDirichletBCMulti.cc \
	TestDirichletBCMultiCases.cc \
	TestDirichletBoundary.cc \
	TestDirichletBoundaryCases.cc \
	TestNeumann.cc \
	TestNeumannCases.cc \
	TestPointForce.cc \
	TestPointForceCases.cc \
	test_bc.cc


noinst_HEADERS = \
	TestBoundaryCondition.hh \
	TestBoundaryConditionPoints.hh \
	TestTimeDependent.hh \
	TestTimeDependentPoints.hh \
	TestBoundaryMesh.hh \
	TestBoundaryMeshCases.hh \
	TestAbsorbingDampers.hh \
	TestAbsorbingDampersCases.hh \
	TestDirichletBC.hh \
	TestDirichletBCCases.hh \
	TestDirichletBCMulti.hh \
	TestDirichletBCMultiCases.hh \
	TestDirichletBoundary.hh \
	TestDirichletBoundaryCases.hh \
	TestNeumann.hh \
	TestNeumannCases.hh \
	TestPointForce.hh \
	TestPointForceCases.hh

# Source files associated with testing data
testbc_SOURCES += \
	data/BoundaryMeshData.cc \
	data/BoundaryMeshDataTri3.cc \
	data/BoundaryMeshDataQuad4.cc \
	data/BoundaryMeshDataTet4.cc \
	data/BoundaryMeshDataHex8.cc \
	data/AbsorbingDampersData.cc \
	data/AbsorbingDampersDataTri3.cc \
	data/AbsorbingDampersDataQuad4.cc \
	data/AbsorbingDampersDataTet4.cc \
	data/AbsorbingDampersDataHex8.cc \
	data/DirichletData.cc \
	data/DirichletDataTri3.cc \
	data/DirichletDataQuad4.cc \
	data/DirichletDataTet4.cc \
	data/DirichletDataHex8.cc \
	data/DirichletDataMulti.cc \
	data/DirichletDataMultiTri3.cc \
	data/DirichletDataMultiTet4.cc \
	data/NeumannData.cc \
	data/NeumannDataTri3.cc \
	data/NeumannDataQuad4.cc \
	data/NeumannDataTet4.cc \
	data/NeumannDataHex8.cc \
	data/PointForceData.cc \
	data/PointForceDataTri3.cc \
	data/PointForceDataQuad4.cc \
	data/PointForceDataTet4.cc \
	data/PointForceDataHex8.cc

noinst_HEADERS += \
	data/BoundaryMeshData.hh \
	data/BoundaryMeshDataTri3.hh \
	data/BoundaryMeshDataQuad4.hh \
	data/BoundaryMeshDataTet4.hh \
	data/BoundaryMeshDataHex8.hh \
	data/AbsorbingDampersData.hh \
	data/AbsorbingDampersDataTri3.hh \
	data/AbsorbingDampersDataQuad4.hh \
	data/AbsorbingDampersDataTet4.hh \
	data/AbsorbingDampersDataHex8.hh \
	data/DirichletData.hh \
	data/DirichletDataTri3.hh \
	data/DirichletDataQuad4.hh \
	data/DirichletDataTet4.hh \
	data/DirichletDataHex8.hh \
	data/DirichletDataMulti.hh \
	data/DirichletDataMultiTri3.hh \
	data/DirichletDataMultiTet4.hh \
	data/NeumannData.hh \
	data/NeumannDataTri3.hh \
	data/NeumannDataQuad4.hh \
	data/NeumannDataTet4.hh \
	data/NeumannDataHex8.hh \
	data/PointForceData.hh \
	data/PointForceDataTri3.hh \
	data/PointForceDataQuad4.hh \
	data/PointForceDataTet4.hh \
	data/PointForceDataHex8.hh


AM_CPPFLAGS += \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES)

testbc_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testbc_LDADD = \
	-lcppunit \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testbc_LDADD += -lnetcdf_c++ -lnetcdf
endif


leakcheck: testbc
	valgrind --log-file=valgrind_bc.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/testbc

# End of file 
