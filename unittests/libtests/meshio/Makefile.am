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

subpackage = meshio
include $(top_srcdir)/subpackage.am

TESTS = testmeshio

check_PROGRAMS = testmeshio

# Primary source files
testmeshio_SOURCES = \
	TestMeshIOAscii.cc \
	test_meshio.cc

noinst_HEADERS = \
	TestMeshIOAscii.hh

# Source files associated with testing data
testmeshio_SOURCES += 

noinst_HEADERS += 

testmeshio_LDFLAGS = $(PETSC_LIB) \
	-ljournal

INCLUDES += $(PETSC_INCLUDE)
AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS)

testmeshio_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/meshio/libpylithmeshio.la

# End of file 
