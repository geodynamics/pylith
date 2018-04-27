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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = materials
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testmaterials

check_PROGRAMS = testmaterials

# Primary source files
testmaterials_SOURCES = \
	TestMaterial.cc \
	TestIsotropicLinearElasticityPlaneStrain.cc \
	TestIsotropicLinearElasticityPlaneStrain_UniformStrain.cc \
	TestIsotropicLinearElasticityPlaneStrain_GravityRefState.cc \
	TestIsotropicLinearMaxwellPlaneStrain.cc \
	TestIsotropicLinearGenMaxwellPlaneStrain.cc \
	TestIsotropicLinearMaxwellPlaneStrain_LinearStrain.cc \
	TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain.cc \
	test_materials.cc


noinst_HEADERS = \
	TestMaterial.hh \
	TestIsotropicLinearElasticityPlaneStrain.hh \
	TestIsotropicLinearMaxwellPlaneStrain.hh \
	TestIsotropicLinearGenMaxwellPlaneStrain.hh

AM_CPPFLAGS += \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES)

testmaterials_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testmaterials_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testmaterials_LDADD += -lnetcdf
endif


# End of file
