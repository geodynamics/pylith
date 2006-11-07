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

lib_LTLIBRARIES = libpylithmeshio.la

libpylithmeshio_la_SOURCES = \
	MeshIO.cc \
	MeshIOAscii.cc

subpkginclude_HEADERS = \
	MeshIO.hh \
	MeshIO.icc \
	MeshIOAscii.hh \
	MeshIOAscii.icc

noinst_HEADERS =

INCLUDES += $(PETSC_INCLUDE)
AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS)

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# version
# $Id$

# End of file 
