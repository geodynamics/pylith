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
	test_bc.cc

noinst_HEADERS = \
	TestBoundaryCondition.hh \
	TestDirichlet.hh \
	TestDirichletLine2.hh \
	TestDirichletLine2b.hh \
	TestDirichletTri3.hh \
	TestDirichletQuad4.hh \
	TestDirichletTet4.hh \
	TestDirichletHex8.hh \
	TestDirichletMulti.hh \
	TestDirichletMultiTri3.hh

# Source files associated with testing data
testbc_SOURCES += \
	data/DirichletData.cc \
	data/DirichletDataLine2.cc \
	data/DirichletDataLine2b.cc \
	data/DirichletDataTri3.cc \
	data/DirichletDataQuad4.cc \
	data/DirichletDataTet4.cc \
	data/DirichletDataHex8.cc \
	data/DirichletDataMulti.cc \
	data/DirichletDataMultiTri3.cc

noinst_HEADERS += \
	data/DirichletData.hh \
	data/DirichletDataLine2.hh \
	data/DirichletDataLine2b.hh \
	data/DirichletDataTri3.hh \
	data/DirichletDataQuad4.hh \
	data/DirichletDataTet4.hh \
	data/DirichletDataHex8.hh \
	data/DirichletDataMulti.hh \
	data/DirichletDataMultiTri3.hh

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
