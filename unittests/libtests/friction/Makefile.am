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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
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
	test_friction.cc \
	TestFrictionModel.cc \
	TestStaticFriction.cc \
	TestSlipWeakening.cc \
	TestRateStateAgeing.cc \
	TestTimeWeakening.cc


noinst_HEADERS = \
	TestFrictionModel.hh \
	TestStaticFriction.hh \
	TestSlipWeakening.hh \
	TestRateStateAgeing.hh \
	TestTimeWeakening.hh

# Source files associated with testing data
testfriction_SOURCES += \
	data/FrictionModelData.cc \
	data/StaticFrictionData.cc \
	data/SlipWeakeningData.cc \
	data/RateStateAgeingData.cc \
	data/TimeWeakeningData.cc


noinst_HEADERS += \
	data/FrictionModelData.hh \
	data/StaticFrictionData.hh \
	data/SlipWeakeningData.hh \
	data/RateStateAgeingData.hh \
	data/TimeWeakeningData.hh

AM_CPPFLAGS = \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES)

testfriction_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testfriction_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testfriction_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
