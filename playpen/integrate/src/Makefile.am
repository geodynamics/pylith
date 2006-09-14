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
MESHIO_INCLUDE = -I../../meshio/src
MESHIO_LIB     = ../../meshio/src/MeshIO.o ../../meshio/src/MeshIOAscii.o

bin_PROGRAMS = testintegrate

testintegrate_SOURCES = \
	elemVector.cc \
	testintegrate.cc

noinst_HEADERS = \
	elemVector.hh

testintegrate_LDADD = $(MESHIO_LIB) $(PETSC_LIB)

INCLUDES = $(PETSC_INCLUDE) $(MESHIO_INCLUDE)

# version
# $Id$

# End of file 
