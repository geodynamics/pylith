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

SUBDIRS = data

TESTS = testmaterials.py
dist_check_SCRIPTS = testmaterials.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestHomogeneous.py \
	TestMaterial.py \
	TestElasticMaterial.py \
	TestElasticStrain1D.py \
	TestElasticStress1D.py \
	TestElasticPlaneStrain.py \
	TestElasticPlaneStress.py \
	TestElasticIsotropic3D.py \
	TestMaxwellIsotropic3D.py \
	TestGenMaxwellIsotropic3D.py


# End of file 
