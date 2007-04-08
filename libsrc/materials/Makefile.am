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
	ElasticStrain1D.hh \
	ElasticStrain1D.icc \
	ElasticStress1D.hh \
	ElasticStress1D.icc \
	ElasticIsotropic3D.hh \
	ElasticIsotropic3D.icc \
	ElasticMaterial.hh \
	ElasticPlaneStrain.hh \
	ElasticPlaneStrain.icc \
	ElasticPlaneStress.hh \
	ElasticPlaneStress.icc \
	Material.hh \
	Material.icc

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
