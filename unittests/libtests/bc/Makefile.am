# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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
	TestAbsorbingDampers.cc \
	TestAbsorbingDampersLine2.cc \
	TestAbsorbingDampersTri3.cc \
	TestAbsorbingDampersQuad4.cc \
	TestAbsorbingDampersTet4.cc \
	TestAbsorbingDampersHex8.cc \
	TestBoundaryCondition.cc \
	TestDirichletPoints.cc \
	TestDirichletPointsLine2.cc \
	TestDirichletPointsLine2b.cc \
	TestDirichletPointsTri3.cc \
	TestDirichletPointsQuad4.cc \
	TestDirichletPointsTet4.cc \
	TestDirichletPointsHex8.cc \
	TestDirichletPointsMulti.cc \
	TestDirichletPointsMultiTri3.cc \
	TestNeumann.cc \
	TestNeumannLine2.cc \
	TestNeumannTri3.cc \
	TestNeumannQuad4.cc \
	TestNeumannTet4.cc \
	TestNeumannHex8.cc \
	test_bc.cc

noinst_HEADERS = \
	TestAbsorbingDampers.hh \
	TestAbsorbingDampersLine2.hh \
	TestAbsorbingDampersTri3.hh \
	TestAbsorbingDampersQuad4.hh \
	TestAbsorbingDampersTet4.hh \
	TestAbsorbingDampersHex8.hh \
	TestBoundaryCondition.hh \
	TestDirichletPoints.hh \
	TestDirichletPointsLine2.hh \
	TestDirichletPointsLine2b.hh \
	TestDirichletPointsTri3.hh \
	TestDirichletPointsQuad4.hh \
	TestDirichletPointsTet4.hh \
	TestDirichletPointsHex8.hh \
	TestDirichletPointsMulti.hh \
	TestDirichletPointsMultiTri3.hh \
	TestNeumann.hh \
	TestNeumannLine2.hh \
	TestNeumannTri3.hh \
	TestNeumannQuad4.hh \
	TestNeumannTet4.hh \
	TestNeumannHex8.hh

# Source files associated with testing data
testbc_SOURCES += \
	data/AbsorbingDampersData.cc \
	data/AbsorbingDampersDataLine2.cc \
	data/AbsorbingDampersDataTri3.cc \
	data/AbsorbingDampersDataQuad4.cc \
	data/AbsorbingDampersDataTet4.cc \
	data/AbsorbingDampersDataHex8.cc \
	data/DirichletPointsData.cc \
	data/DirichletPointsDataLine2.cc \
	data/DirichletPointsDataLine2b.cc \
	data/DirichletPointsDataTri3.cc \
	data/DirichletPointsDataQuad4.cc \
	data/DirichletPointsDataTet4.cc \
	data/DirichletPointsDataHex8.cc \
	data/DirichletPointsDataMulti.cc \
	data/DirichletPointsDataMultiTri3.cc \
	data/NeumannData.cc \
	data/NeumannDataLine2.cc \
	data/NeumannDataTri3.cc \
	data/NeumannDataQuad4.cc \
	data/NeumannDataTet4.cc \
	data/NeumannDataHex8.cc

noinst_HEADERS += \
	data/AbsorbingDampersData.hh \
	data/AbsorbingDampersDataLine2.hh \
	data/AbsorbingDampersDataTri3.hh \
	data/AbsorbingDampersDataQuad4.hh \
	data/AbsorbingDampersDataTet4.hh \
	data/AbsorbingDampersDataHex8.hh \
	data/DirichletPointsData.hh \
	data/DirichletPointsDataLine2.hh \
	data/DirichletPointsDataLine2b.hh \
	data/DirichletPointsDataTri3.hh \
	data/DirichletPointsDataQuad4.hh \
	data/DirichletPointsDataTet4.hh \
	data/DirichletPointsDataHex8.hh \
	data/DirichletPointsDataMulti.hh \
	data/DirichletPointsDataMultiTri3.hh \
	data/NeumannData.hh \
	data/NeumannDataLine2.hh \
	data/NeumannDataTri3.hh \
	data/NeumannDataQuad4.hh \
	data/NeumannDataTet4.hh \
	data/NeumannDataHex8.hh

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testbc_LDADD = \
	-lcppunit \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testbc_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
