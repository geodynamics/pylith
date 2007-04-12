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

subpackage = materials
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testmaterials

check_PROGRAMS = testmaterials

# Primary source files
testmaterials_SOURCES = \
	TestMaterial.cc \
	TestElasticMaterial.cc \
	TestElasticIsotropic3D.cc \
	TestElasticPlaneStrain.cc \
	TestElasticPlaneStress.cc \
	TestElasticStrain1D.cc \
	TestElasticStress1D.cc \
	test_materials.cc

noinst_HEADERS = \
	TestMaterial.hh

# Source files associated with testing data
testmaterials_SOURCES += \
	data/MaterialData.cc \
	data/ElasticMaterialData.cc \
	data/ElasticIsotropic3DData.cc \
	data/ElasticPlaneStrainData.cc \
	data/ElasticPlaneStressData.cc \
	data/ElasticStrain1DData.cc \
	data/ElasticStress1DData.cc

noinst_HEADERS += \
	data/MaterialData.hh \
	data/ElasticMaterialData.hh \
	data/ElasticIsotropic3DData.hh \
	data/ElasticPlaneStrainData.hh \
	data/ElasticPlaneStressData.hh \
	data/ElasticStrain1DData.hh \
	data/ElasticStress1DData.hh \
	data/header.hh

testmaterials_LDFLAGS = $(PETSC_LIB)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testmaterials_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata


# End of file 
