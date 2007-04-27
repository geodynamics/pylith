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
	TestFault.cc \
	TestFaultCohesive.cc \
	TestFaultCohesiveKin.cc \
	test_faults.cc

noinst_HEADERS = \
	TestFault.hh

# Source files associated with testing data
testfaults_SOURCES += \
	data/CohesiveData.cc \
	data/CohesiveDataLine2.cc \
	data/CohesiveLagrangeDataLine2.cc \
	data/CohesiveDataTri3.cc \
	data/CohesiveLagrangeDataTri3.cc \
	data/CohesiveDataQuad4.cc \
	data/CohesiveLagrangeDataQuad4.cc \
	data/CohesiveDataHex8.cc \
	data/CohesiveLagrangeDataHex8.cc \
	data/CohesiveLagrangeDataTet4.cc \
	data/CohesiveDataTet4.cc

noinst_HEADERS += \
	data/CohesiveData.hh \
	data/CohesiveDataLine2.hh \
	data/CohesiveLagrangeDataLine2.hh \
	data/CohesiveDataTri3.hh \
	data/CohesiveLagrangeDataTri3.hh \
	data/CohesiveDataQuad4.hh \
	data/CohesiveLagrangeDataQuad4.hh \
	data/CohesiveDataHex8.hh \
	data/CohesiveLagrangeDataHex8.hh \
	data/CohesiveDataTet4.hh \
	data/CohesiveLagrangeDataTet4.hh

testfaults_LDFLAGS = $(PETSC_LIB)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testfaults_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata


# End of file 
