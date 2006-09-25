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
	Quadrature.cc \
	Quadrature1D.cc \
	Quadrature1Din2D.cc \
	Quadrature1Din3D.cc \
	Quadrature2D.cc

subpkginclude_HEADERS = \
	Quadrature.hh \
	Quadrature.icc \
	Quadrature1D.hh \
	Quadrature1D.icc \
	Quadrature1Din2D.hh \
	Quadrature1Din2D.icc \
	Quadrature1Din3D.hh \
	Quadrature1Din3D.icc \
	Quadrature2D.hh \
	Quadrature2D.icc

noinst_HEADERS =

INCLUDES += $(PETSC_INCLUDE)

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# version
# $Id$

# End of file 
