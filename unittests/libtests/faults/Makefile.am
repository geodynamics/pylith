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

subpackage = faults
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testfaults

check_PROGRAMS = testfaults

# Primary source files
testfaults_SOURCES = \
	TestBruneSlipFn.cc \
	TestEqKinSrc.cc \
	TestFault.cc \
	TestFaultCohesive.cc \
	TestFaultCohesiveKin.cc \
	TestFaultCohesiveKinLine2.cc \
	TestFaultCohesiveKinTri3.cc \
	TestFaultCohesiveKinTri3d.cc \
	TestFaultCohesiveKinQuad4.cc \
	TestFaultCohesiveKinQuad4e.cc \
	TestFaultCohesiveKinTet4.cc \
	TestFaultCohesiveKinHex8.cc \
	test_faults.cc

noinst_HEADERS = \
	TestBruneSlipFn.hh \
	TestEqKinSrc.hh \
	TestFault.hh \
	TestFaultCohesive.hh \
	TestFaultCohesiveKin.hh \
	TestFaultCohesiveKinLine2.hh \
	TestFaultCohesiveKinTri3.hh \
	TestFaultCohesiveKinTri3d.hh \
	TestFaultCohesiveKinQuad4.hh \
	TestFaultCohesiveKinQuad4e.hh \
	TestFaultCohesiveKinTet4.hh \
	TestFaultCohesiveKinHex8.hh

# Source files associated with testing data
testfaults_SOURCES += \
	data/CohesiveData.cc \
	data/CohesiveDataLine2.cc \
	data/CohesiveDataLine2Lagrange.cc \
	data/CohesiveDataTri3.cc \
	data/CohesiveDataTri3b.cc \
	data/CohesiveDataTri3c.cc \
	data/CohesiveDataTri3d.cc \
	data/CohesiveDataTri3Lagrange.cc \
	data/CohesiveDataQuad4.cc \
	data/CohesiveDataQuad4b.cc \
	data/CohesiveDataQuad4c.cc \
	data/CohesiveDataQuad4d.cc \
	data/CohesiveDataQuad4e.cc \
	data/CohesiveDataQuad4Lagrange.cc \
	data/CohesiveDataTet4.cc \
	data/CohesiveDataTet4b.cc \
	data/CohesiveDataTet4c.cc \
	data/CohesiveDataTet4d.cc \
	data/CohesiveDataTet4Lagrange.cc \
	data/CohesiveDataHex8.cc \
	data/CohesiveDataHex8b.cc \
	data/CohesiveDataHex8c.cc \
	data/CohesiveDataHex8d.cc \
	data/CohesiveDataHex8e.cc \
	data/CohesiveDataHex8f.cc \
	data/CohesiveDataHex8Lagrange.cc \
	data/CohesiveKinData.cc \
	data/CohesiveKinDataLine2.cc \
	data/CohesiveKinDataTri3.cc \
	data/CohesiveKinDataTri3d.cc \
	data/CohesiveKinDataQuad4.cc \
	data/CohesiveKinDataQuad4e.cc \
	data/CohesiveKinDataTet4.cc \
	data/CohesiveKinDataHex8.cc

noinst_HEADERS += \
	data/CohesiveData.hh \
	data/CohesiveDataLine2.hh \
	data/CohesiveDataLine2Lagrange.hh \
	data/CohesiveDataTri3.hh \
	data/CohesiveDataTri3b.hh \
	data/CohesiveDataTri3c.hh \
	data/CohesiveDataTri3d.hh \
	data/CohesiveDataTri3Lagrange.hh \
	data/CohesiveDataQuad4.hh \
	data/CohesiveDataQuad4b.hh \
	data/CohesiveDataQuad4c.hh \
	data/CohesiveDataQuad4d.hh \
	data/CohesiveDataQuad4e.hh \
	data/CohesiveDataQuad4Lagrange.hh \
	data/CohesiveDataTet4.hh \
	data/CohesiveDataTet4b.hh \
	data/CohesiveDataTet4c.hh \
	data/CohesiveDataTet4d.hh \
	data/CohesiveDataTet4Lagrange.hh \
	data/CohesiveDataHex8.hh \
	data/CohesiveDataHex8b.hh \
	data/CohesiveDataHex8c.hh \
	data/CohesiveDataHex8d.hh \
	data/CohesiveDataHex8e.hh \
	data/CohesiveDataHex8f.hh \
	data/CohesiveDataHex8Lagrange.hh \
	data/CohesiveKinData.hh \
	data/CohesiveKinDataLine2.hh \
	data/CohesiveKinDataTri3.hh \
	data/CohesiveKinDataTri3d.hh \
	data/CohesiveKinDataQuad4.hh \
	data/CohesiveKinDataQuad4e.hh \
	data/CohesiveKinDataTet4.hh \
	data/CohesiveKinDataHex8.hh

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testfaults_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testfaults_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
