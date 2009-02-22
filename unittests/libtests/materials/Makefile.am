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
	TestMetadata.cc \
	TestMaterial.cc \
	TestElasticMaterial.cc \
	TestElasticIsotropic3D.cc \
	test_materials.cc

#	TestElasticPlaneStrain.cc \
#	TestElasticPlaneStress.cc \
#	TestElasticStrain1D.cc \
#	TestElasticStress1D.cc \
#	TestGenMaxwellIsotropic3D.cc \
#	TestMaxwellIsotropic3D.cc 



noinst_HEADERS = \
	TestMaterial.hh \
	TestElasticMaterial.hh \
	TestElasticStrain1D.hh \
	TestElasticStress1D.hh \
	TestElasticPlaneStrain.hh \
	TestElasticPlaneStress.hh \
	TestElasticIsotropic3D.hh \
	TestGenMaxwellIsotropic3D.hh \
	TestMaxwellIsotropic3D.hh

# Source files associated with testing data
testmaterials_SOURCES += \
	data/MaterialData.cc \
	data/ElasticMaterialData.cc \
	data/ElasticIsotropic3DData.cc 

#	data/ElasticPlaneStrainData.cc \
#	data/ElasticPlaneStressData.cc \
#	data/ElasticStrain1DData.cc \
#	data/ElasticStress1DData.cc \
#	data/GenMaxwellIsotropic3DElasticData.cc \
#	data/GenMaxwellIsotropic3DTimeDepData.cc \
#	data/MaxwellIsotropic3DElasticData.cc \
#	data/MaxwellIsotropic3DTimeDepData.cc

noinst_HEADERS += \
	data/MaterialData.hh \
	data/ElasticMaterialData.hh \
	data/ElasticIsotropic3DData.hh \
	data/ElasticPlaneStrainData.hh \
	data/ElasticPlaneStressData.hh \
	data/ElasticStrain1DData.hh \
	data/ElasticStress1DData.hh \
	data/GenMaxwellIsotropic3DElasticData.hh \
	data/GenMaxwellIsotropic3DTimeDepData.hh \
	data/MaxwellIsotropic3DElasticData.hh \
	data/MaxwellIsotropic3DTimeDepData.hh \
	data/header.hh

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testmaterials_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testmaterials_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
