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
	MeshIO.cc \
	MeshIOAscii.cc \
	testascii.cc

#	PetscMesh.cc
#	HDF5.cc
#	MeshIOHDF5.cc

noinst_HEADERS = \
	MeshIO.hh \
	MeshIOAscii.hh

#	PetscMesh.hh
#	MeshIOHDF5.hh \
#	HDF5.hh

#HDF5_LIBS = -lhdf5
testascii_LDADD = $(PETSC_LIB) $(HDF5_LIBS)

INCLUDES = $(PETSC_INCLUDE)

# version
# $Id$

# End of file 
