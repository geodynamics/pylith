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

subpackage = materials
include $(top_srcdir)/subpackage.am

lib_LTLIBRARIES = libpylithmaterials.la

libpylithmaterials_la_SOURCES = \
	ElasticIsotropic3D.cc \
	ElasticMaterial.cc

subpkginclude_HEADERS = \
	ElasticIsotropic3D.hh \
	ElasticIsotropic3D.icc \
	ElasticMaterial.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# version
# $Id$

# End of file 
