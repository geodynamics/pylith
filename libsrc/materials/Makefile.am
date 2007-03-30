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
	ElasticIsotropic1D.hh \
	ElasticIsotropic1D.icc \
	ElasticIsotropic3D.hh \
	ElasticIsotropic3D.icc \
	ElasticMaterial.hh \
	ElasticMaterial.icc \
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
