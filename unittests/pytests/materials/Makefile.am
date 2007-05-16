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

noinst_PYTHON = \
	TestHomogeneous.py \
	TestMaterial.py \
	TestElasticIsotropic3D.py

TESTS_ENVIRONMENT = $(PYTHON)


# End of file 
