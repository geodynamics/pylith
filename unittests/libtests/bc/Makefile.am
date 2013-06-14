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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = bc
include $(top_srcdir)/subpackage.am

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
	TestBoundaryMeshTri3.cc \
	TestBoundaryMeshQuad4.cc \
	TestBoundaryMeshTet4.cc \
	TestBoundaryMeshHex8.cc \
	TestAbsorbingDampers.cc \
	TestAbsorbingDampersLine2.cc \
	TestAbsorbingDampersTri3.cc \
	TestAbsorbingDampersQuad4.cc \
	TestAbsorbingDampersTet4.cc \
	TestAbsorbingDampersHex8.cc \
	TestDirichletBC.cc \
	TestDirichletBCLine2.cc \
	TestDirichletBCLine2b.cc \
	TestDirichletBCTri3.cc \
	TestDirichletBCQuad4.cc \
	TestDirichletBCTet4.cc \
	TestDirichletBCHex8.cc \
	TestDirichletBCMulti.cc \
	TestDirichletBCMultiTri3.cc \
	TestDirichletBCMultiTet4.cc \
	TestDirichletBoundary.cc \
	TestDirichletBoundaryTri3.cc \
	TestDirichletBoundaryQuad4.cc \
	TestDirichletBoundaryTet4.cc \
	TestDirichletBoundaryHex8.cc \
	TestNeumann.cc \
	TestNeumannLine2.cc \
	TestNeumannTri3.cc \
	TestNeumannQuad4.cc \
	TestNeumannTet4.cc \
	TestNeumannHex8.cc \
	TestPointForce.cc \
	TestPointForceLine2.cc \
	TestPointForceTri3.cc \
	TestPointForceQuad4.cc \
	TestPointForceTet4.cc \
	TestPointForceHex8.cc \
	test_bc.cc


noinst_HEADERS = \
	TestBoundaryCondition.hh \
	TestBoundaryConditionPoints.hh \
	TestTimeDependent.hh \
	TestTimeDependentPoints.hh \
	TestBoundaryMesh.hh \
	TestBoundaryMeshTri3.hh \
	TestBoundaryMeshQuad4.hh \
	TestBoundaryMeshTet4.hh \
	TestBoundaryMeshHex8.hh \
	TestAbsorbingDampers.hh \
	TestAbsorbingDampersLine2.hh \
	TestAbsorbingDampersTri3.hh \
	TestAbsorbingDampersQuad4.hh \
	TestAbsorbingDampersTet4.hh \
	TestAbsorbingDampersHex8.hh \
	TestDirichletBC.hh \
	TestDirichletBCLine2.hh \
	TestDirichletBCLine2b.hh \
	TestDirichletBCTri3.hh \
	TestDirichletBCQuad4.hh \
	TestDirichletBCTet4.hh \
	TestDirichletBCHex8.hh \
	TestDirichletBCMulti.hh \
	TestDirichletBCMultiTri3.hh \
	TestDirichletBCMultiTet4.hh \
	TestDirichletBoundary.hh \
	TestDirichletBoundaryTri3.hh \
	TestDirichletBoundaryQuad4.hh \
	TestDirichletBoundaryTet4.hh \
	TestDirichletBoundaryHex8.hh \
	TestNeumann.hh \
	TestNeumannLine2.hh \
	TestNeumannTri3.hh \
	TestNeumannQuad4.hh \
	TestNeumannTet4.hh \
	TestNeumannHex8.hh \
	TestPointForce.hh \
	TestPointForceLine2.hh \
	TestPointForceTri3.hh \
	TestPointForceQuad4.hh \
	TestPointForceTet4.hh \
	TestPointForceHex8.hh

# Source files associated with testing data
testbc_SOURCES += \
	data/BoundaryMeshData.cc \
	data/BoundaryMeshDataTri3.cc \
	data/BoundaryMeshDataQuad4.cc \
	data/BoundaryMeshDataTet4.cc \
	data/BoundaryMeshDataHex8.cc \
	data/AbsorbingDampersData.cc \
	data/AbsorbingDampersDataLine2.cc \
	data/AbsorbingDampersDataTri3.cc \
	data/AbsorbingDampersDataQuad4.cc \
	data/AbsorbingDampersDataTet4.cc \
	data/AbsorbingDampersDataHex8.cc \
	data/DirichletData.cc \
	data/DirichletDataLine2.cc \
	data/DirichletDataLine2b.cc \
	data/DirichletDataTri3.cc \
	data/DirichletDataQuad4.cc \
	data/DirichletDataTet4.cc \
	data/DirichletDataHex8.cc \
	data/DirichletDataMulti.cc \
	data/DirichletDataMultiTri3.cc \
	data/DirichletDataMultiTet4.cc \
	data/NeumannData.cc \
	data/NeumannDataLine2.cc \
	data/NeumannDataTri3.cc \
	data/NeumannDataQuad4.cc \
	data/NeumannDataTet4.cc \
	data/NeumannDataHex8.cc \
	data/PointForceData.cc \
	data/PointForceDataLine2.cc \
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
	data/AbsorbingDampersDataLine2.hh \
	data/AbsorbingDampersDataTri3.hh \
	data/AbsorbingDampersDataQuad4.hh \
	data/AbsorbingDampersDataTet4.hh \
	data/AbsorbingDampersDataHex8.hh \
	data/DirichletData.hh \
	data/DirichletDataLine2.hh \
	data/DirichletDataLine2b.hh \
	data/DirichletDataTri3.hh \
	data/DirichletDataQuad4.hh \
	data/DirichletDataTet4.hh \
	data/DirichletDataHex8.hh \
	data/DirichletDataMulti.hh \
	data/DirichletDataMultiTri3.hh \
	data/DirichletDataMultiTet4.hh \
	data/NeumannData.hh \
	data/NeumannDataLine2.hh \
	data/NeumannDataTri3.hh \
	data/NeumannDataQuad4.hh \
	data/NeumannDataTet4.hh \
	data/NeumannDataHex8.hh \
	data/PointForceData.hh \
	data/PointForceDataLine2.hh \
	data/PointForceDataTri3.hh \
	data/PointForceDataQuad4.hh \
	data/PointForceDataTet4.hh \
	data/PointForceDataHex8.hh


AM_CPPFLAGS = \
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


# End of file 
