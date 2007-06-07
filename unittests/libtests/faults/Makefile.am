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
	test_faults.cc

noinst_HEADERS = \
	TestBruneSlipFn.hh \
	TestEqKinSrc.hh \
	TestFault.hh \
	TestFaultCohesive.hh \
	TestFaultCohesiveKin.hh

# Source files associated with testing data
testfaults_SOURCES += \
	data/CohesiveData.cc \
	data/CohesiveDataLine2.cc \
	data/CohesiveDataLine2Lagrange.cc \
	data/CohesiveDataTri3.cc \
	data/CohesiveDataTri3Lagrange.cc \
	data/CohesiveDataQuad4.cc \
	data/CohesiveDataQuad4Lagrange.cc \
	data/CohesiveDataHex8.cc \
	data/CohesiveDataHex8Lagrange.cc \
	data/CohesiveDataTet4.cc \
	data/CohesiveDataTet4Lagrange.cc

noinst_HEADERS += \
	data/CohesiveData.hh \
	data/CohesiveDataLine2.hh \
	data/CohesiveDataLine2Lagrange.hh \
	data/CohesiveDataTri3.hh \
	data/CohesiveDataTri3Lagrange.hh \
	data/CohesiveDataQuad4.hh \
	data/CohesiveDataQuad4Lagrange.hh \
	data/CohesiveDataHex8.hh \
	data/CohesiveDataHex8Lagrange.hh \
	data/CohesiveDataTet4.hh \
	data/CohesiveDataTet4Lagrange.hh

testfaults_LDFLAGS = $(PETSC_LIB) $(PYTHON_BLDLIBRARY)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testfaults_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata

if ENABLE_CUBIT
  testfaults_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
