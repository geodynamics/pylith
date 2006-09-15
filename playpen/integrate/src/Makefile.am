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
CIG_INCLUDE = -I/PETSc3/geoframe/packages/pylith-new/include
CIG_LIB     = -L/PETSc3/geoframe/packages/pylith-new/lib -lpylithmeshio

bin_PROGRAMS = testintegrate

testintegrate_SOURCES = \
	elemVector.cc \
	testintegrate.cc

noinst_HEADERS = \
	elemVector.hh

testintegrate_LDADD = $(CIG_LIB) $(PETSC_LIB)

INCLUDES = $(CIG_INCLUDE) $(PETSC_INCLUDE)

# version
# $Id$

# End of file 
