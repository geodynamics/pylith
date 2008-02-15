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
check_SCRIPTS = testmaterials.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestBiMaterial.py \
	TestElasticIsotropic3D.py \
	TestElasticPlaneStrain.py \
	TestElasticPlaneStress.py \
	TestElasticStrain1D.py \
	TestElasticStress1D.py \
	TestGenMaxwellIsotropic3D.py \
	TestHomogeneous.py \
	TestMaterial.py \
	TestMaxwellIsotropic3D.py


# End of file 
