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

TESTS = testmaterials

check_PROGRAMS = testmaterials

# Primary source files
testmaterials_SOURCES = \
	TestMaterial.cc \
	test_materials.cc

noinst_HEADERS = \
	TestMaterial.hh

# Source files associated with testing data
testmaterials_SOURCES += \
	data/MaterialData.cc \
	data/ElasticIsotropic3DData.cc

noinst_HEADERS += \
	data/MaterialData.hh \
	data/ElasticIsotropic3DData.hh

testmaterials_LDFLAGS = $(PETSC_LIB)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testmaterials_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata

noinst_tmp =

CLEANFILES = $(noinst_tmp)

# End of file 
