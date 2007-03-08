# -*- Makefile -*-
#
CIG_INCLUDE = -I/PETSc3/cig/pylith3d/include
CIG_LIB     = -L/PETSc3/cig/pylith3d/lib -lpylith

bin_PROGRAMS = testcohesive

testcohesive_SOURCES = \
	testcohesive.cc

#noinst_HEADERS = \
#	elemVector.hh

testcohesive_LDADD = $(CIG_LIB) $(PETSC_LIB)

INCLUDES = $(CIG_INCLUDE) $(PETSC_INCLUDE)

# version
# $Id$

# End of file 
