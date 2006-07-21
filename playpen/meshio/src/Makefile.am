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
	MeshIOAscii.cc \
	PetscMesh.cc \
	testascii.cc

noinst_HEADERS = \
	MeshIOAscii.hh \
	PetscMesh.hh

testascii_LDADD = $(PETSC_LIB)

INCLUDES = $(PETSC_INCLUDE)

# version
# $Id$

# End of file 
