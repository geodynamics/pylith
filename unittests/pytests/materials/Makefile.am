# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
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
	TestMaxwellPlaneStrain.py \
	TestGenMaxwellIsotropic3D.py \
	TestGenMaxwellPlaneStrain.py \
	TestGenMaxwellQpQsIsotropic3D.py \
	TestPowerLaw3D.py \
	TestDruckerPrager3D.py


# End of file 
