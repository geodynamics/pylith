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

bin_PROGRAMS = testascii

testascii_SOURCES = \
	HDF5.cc

#	MeshIO.cc \
#	MeshIOAscii.cc \
#	MeshIOHDF5.cc \
#	PetscMesh.cc \
#	testascii.cc

noinst_HEADERS = \
	HDF5.hh \
	MeshIO.hh \
	MeshIOAscii.hh \
	MeshIOHDF5.hh \
	PetscMesh.hh

HDF5_LIBS = -lhdf5
testascii_LDADD = $(PETSC_LIB) $(HDF5_LIBS)

INCLUDES = $(PETSC_INCLUDE)

# version
# $Id$

# End of file 
