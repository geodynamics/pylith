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

subpackage = topology
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testtopology

check_PROGRAMS = testtopology

# Primary source files
testtopology_SOURCES = \
	TestMesh.cc \
	TestSubMesh.cc \
	TestMeshOps.cc \
	TestFieldBase.cc \
	TestField.cc \
	TestFieldSubMesh.cc \
	test_topology.cc

#	TestFieldsManager.cc

noinst_HEADERS = \
	TestMesh.hh \
	TestSubMesh.hh \
	TestMeshOps.hh \
	TestFieldBase.hh \
	TestField.hh \
	TestFieldSubMesh.hh \
	TestFieldsManager.hh

# Source files associated with testing data
testtopology_SOURCES += 

noinst_HEADERS += 

AM_CPPFLAGS = \
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testtopology_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testtopology_LDADD += -lnetcdf_c++ -lnetcdf
endif

noinst_tmp =

CLEANFILES = $(noinst_tmp)

# End of file 
