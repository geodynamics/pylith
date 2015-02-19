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
# Copyright (c) 2010-2015 University of California, Davis
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
	TestMetadata.cc \
	TestMaterial.cc \
	TestElasticMaterial.cc \
	TestElasticPlaneStrain.cc \
	TestElasticPlaneStress.cc \
	TestElasticIsotropic3D.cc \
	TestMaxwellIsotropic3D.cc \
	TestMaxwellPlaneStrain.cc \
	TestGenMaxwellIsotropic3D.cc \
	TestGenMaxwellPlaneStrain.cc \
	TestGenMaxwellQpQsIsotropic3D.cc \
	TestPowerLaw3D.cc \
	TestPowerLawPlaneStrain.cc \
	TestDruckerPrager3D.cc \
	TestDruckerPragerPlaneStrain.cc \
	TestEffectiveStress.cc \
	test_materials.cc


noinst_HEADERS = \
	TestMetadata.hh \
	TestMaterial.hh \
	TestElasticMaterial.hh \
	TestElasticPlaneStrain.hh \
	TestElasticPlaneStress.hh \
	TestElasticIsotropic3D.hh \
	TestGenMaxwellIsotropic3D.hh \
	TestGenMaxwellPlaneStrain.hh \
	TestGenMaxwellQpQsIsotropic3D.hh \
	TestMaxwellIsotropic3D.hh \
	TestMaxwellPlaneStrain.hh \
	TestPowerLaw3D.hh \
	TestPowerLawPlaneStrain.hh \
	TestDruckerPrager3D.hh \
	TestDruckerPragerPlaneStrain.hh \
	TestEffectiveStress.hh

# Source files associated with testing data
testmaterials_SOURCES += \
	data/MaterialData.cc \
	data/ElasticMaterialData.cc \
	data/ElasticPlaneStrainData.cc \
	data/ElasticPlaneStressData.cc \
	data/ElasticIsotropic3DData.cc \
	data/MaxwellIsotropic3DElasticData.cc \
	data/MaxwellIsotropic3DTimeDepData.cc \
	data/MaxwellPlaneStrainElasticData.cc \
	data/MaxwellPlaneStrainTimeDepData.cc \
	data/GenMaxwellIsotropic3DElasticData.cc \
	data/GenMaxwellIsotropic3DTimeDepData.cc \
	data/GenMaxwellPlaneStrainElasticData.cc \
	data/GenMaxwellPlaneStrainTimeDepData.cc \
	data/GenMaxwellQpQsIsotropic3DElasticData.cc \
	data/GenMaxwellQpQsIsotropic3DTimeDepData.cc \
	data/PowerLaw3DElasticData.cc \
	data/PowerLaw3DTimeDepData.cc \
	data/PowerLawPlaneStrainElasticData.cc \
	data/PowerLawPlaneStrainTimeDepData.cc \
	data/DruckerPrager3DElasticData.cc \
	data/DruckerPrager3DTimeDepData.cc \
	data/DruckerPragerPlaneStrainElasticData.cc \
	data/DruckerPragerPlaneStrainTimeDepData.cc


noinst_HEADERS += \
	data/MaterialData.hh \
	data/ElasticMaterialData.hh \
	data/ElasticIsotropic3DData.hh \
	data/ElasticPlaneStrainData.hh \
	data/ElasticPlaneStressData.hh \
	data/GenMaxwellIsotropic3DElasticData.hh \
	data/GenMaxwellIsotropic3DTimeDepData.hh \
	data/GenMaxwellPlaneStrainElasticData.hh \
	data/GenMaxwellPlaneStrainTimeDepData.hh \
	data/GenMaxwellQpQsIsotropic3DElasticData.hh \
	data/GenMaxwellQpQsIsotropic3DTimeDepData.hh \
	data/MaxwellIsotropic3DElasticData.hh \
	data/MaxwellIsotropic3DTimeDepData.hh \
	data/MaxwellPlaneStrainElasticData.hh \
	data/MaxwellPlaneStrainTimeDepData.hh \
	data/PowerLaw3DElasticData.hh \
	data/PowerLaw3DTimeDepData.hh \
	data/PowerLawPlaneStrainElasticData.hh \
	data/PowerLawPlaneStrainTimeDepData.hh \
	data/DruckerPrager3DElasticData.hh \
	data/DruckerPrager3DTimeDepData.hh \
	data/DruckerPragerPlaneStrainElasticData.hh \
	data/DruckerPragerPlaneStrainTimeDepData.hh \
	data/header.hh

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
  testmaterials_LDADD += -lnetcdf_c++ -lnetcdf
endif


leakcheck: testmaterials
	valgrind --log-file=valgrind_materials.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/testmaterials


# End of file 
