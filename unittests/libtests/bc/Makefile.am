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
	test_bc.cc

noinst_HEADERS = \
	TestBoundaryCondition.hh \
	TestDirichletLine2.hh \
	TestDirichlet.hh

# Source files associated with testing data
testbc_SOURCES += \
	data/DirichletData.cc \
	data/DirichletDataLine2.cc

noinst_HEADERS += \
	data/DirichletData.hh \
	data/DirichletDataLine2.hh

testbc_LDFLAGS = $(PETSC_LIB) $(PYTHON_BLDLIBRARY)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testbc_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata

if ENABLE_CUBIT
  testbc_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
