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
	TestDirichlet.cc \
	TestDirichletLine2.cc \
	TestDirichletLine2b.cc \
	TestDirichletTri3.cc \
	TestDirichletQuad4.cc \
	TestDirichletTet4.cc \
	TestDirichletHex8.cc \
	TestDirichletMulti.cc \
	TestDirichletMultiTri3.cc \
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
	TestDirichlet.hh \
	TestDirichletLine2.hh \
	TestDirichletLine2b.hh \
	TestDirichletTri3.hh \
	TestDirichletQuad4.hh \
	TestDirichletTet4.hh \
	TestDirichletHex8.hh \
	TestDirichletMulti.hh \
	TestDirichletMultiTri3.hh \
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
	data/DirichletData.cc \
	data/DirichletDataLine2.cc \
	data/DirichletDataLine2b.cc \
	data/DirichletDataTri3.cc \
	data/DirichletDataQuad4.cc \
	data/DirichletDataTet4.cc \
	data/DirichletDataHex8.cc \
	data/DirichletDataMulti.cc \
	data/DirichletDataMultiTri3.cc \
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
	data/DirichletData.hh \
	data/DirichletDataLine2.hh \
	data/DirichletDataLine2b.hh \
	data/DirichletDataTri3.hh \
	data/DirichletDataQuad4.hh \
	data/DirichletDataTet4.hh \
	data/DirichletDataHex8.hh \
	data/DirichletDataMulti.hh \
	data/DirichletDataMultiTri3.hh \
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
