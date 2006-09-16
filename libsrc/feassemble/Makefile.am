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

subpackage = feassemble
include $(top_srcdir)/subpackage.am

lib_LTLIBRARIES = libpylithfeassemble.la

libpylithfeassemble_la_SOURCES = \
	Quadrature.cc

subpkginclude_HEADERS = \
	Quadrature.hh

noinst_HEADERS =

INCLUDES += $(PETSC_INCLUDE)

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# version
# $Id$

# End of file 
