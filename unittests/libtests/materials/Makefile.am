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

TESTS = \
	test_isotropiclinearelasticityplanestrain \
	test_isotropiclinearmaxwellplanestrain \
	test_isotropiclineargenmaxwellplanestrain \
	test_isotropiclinearelasticity3d \
	test_isotropiclinearmaxwell3d \
	test_isotropiclineargenmaxwell3d

check_PROGRAMS = $(TESTS)


test_isotropiclinearelasticityplanestrain_SOURCES = \
	TestMaterial.cc \
	TestIsotropicLinearElasticityPlaneStrain.cc \
	TestIsotropicLinearElasticityPlaneStrain_UniformStrain.cc \
	TestIsotropicLinearElasticityPlaneStrain_GravityRefState.cc \
	test_materials.cc

test_isotropiclinearmaxwellplanestrain_SOURCES = \
	TestMaterial.cc \
	TestIsotropicLinearMaxwellPlaneStrain.cc \
	TestIsotropicLinearMaxwellPlaneStrain_LinearStrain.cc \
	test_materials.cc

test_isotropiclineargenmaxwellplanestrain_SOURCES = \
	TestMaterial.cc \
	TestIsotropicLinearGenMaxwellPlaneStrain.cc \
	TestIsotropicLinearGenMaxwellPlaneStrain_LinearStrain.cc \
	test_materials.cc

test_isotropiclinearelasticity3d_SOURCES = \
	TestMaterial.cc \
	TestIsotropicLinearElasticity3D.cc \
	TestIsotropicLinearElasticity3D_UniformStrain.cc \
	test_materials.cc

test_isotropiclinearmaxwell3d_SOURCES = \
	TestMaterial.cc \
	TestIsotropicLinearMaxwell3D.cc \
	TestIsotropicLinearMaxwell3D_LinearStrain.cc \
	test_materials.cc

test_isotropiclineargenmaxwell3d_SOURCES = \
	TestMaterial.cc \
	TestIsotropicLinearGenMaxwell3D.cc \
	TestIsotropicLinearGenMaxwell3D_LinearStrain.cc \
	test_materials.cc

noinst_HEADERS = \
	TestMaterial.hh \
	TestIsotropicLinearElasticityPlaneStrain.hh \
	TestIsotropicLinearElasticity3D.hh \
	TestIsotropicLinearMaxwellPlaneStrain.hh \
	TestIsotropicLinearGenMaxwellPlaneStrain.hh \
	TestIsotropicLinearMaxwell3D.hh \
	TestIsotropicLinearGenMaxwell3D.hh

AM_CPPFLAGS += \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES)

AM_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  LDADD += -lnetcdf
endif


# End of file
