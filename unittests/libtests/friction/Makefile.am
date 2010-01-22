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

subpackage = friction
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testfriction

check_PROGRAMS = testfriction

# Primary source files
testfriction_SOURCES = \
	test_friction.cc

#	TestFrictionModel.cc \
#	TestStaticFriction.cc


noinst_HEADERS = 

#	TestFrictionModel.hh \
#	TestStaticFriction.hh

# Source files associated with testing data
testfriction_SOURCES += 

#	data/FrictionModelData.cc \
#	data/StaticFrictionData.cc


noinst_HEADERS += 

#	data/FrictionModelData.hh \
#	data/StaticFrictionData.hh

AM_CPPFLAGS = \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testfriction_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testfriction_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testfriction_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
