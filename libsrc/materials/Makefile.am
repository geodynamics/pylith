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

subpkginclude_HEADERS = \
	ElasticIsotropic3D.hh \
	ElasticIsotropic3D.icc \
	ElasticMaterial.hh \
	ElasticMaterial.icc \
	Material.hh \
	Material.icc

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
