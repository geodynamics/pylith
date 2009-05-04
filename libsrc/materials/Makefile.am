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
	ElasticStress1D.hh \
	ElasticIsotropic3D.hh \
	ElasticMaterial.hh \
	ElasticMaterial.icc \
	ElasticPlaneStrain.hh \
	ElasticPlaneStress.hh \
	GenMaxwellIsotropic3D.hh \
	GenMaxwellIsotropic3D.icc \
	MaxwellIsotropic3D.hh \
	MaxwellIsotropic3D.icc \
	Metadata.hh \
	Metadata.icc \
	Material.hh \
	Material.icc \
	ViscoelasticMaxwell.hh \
	materialsfwd.hh


noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
