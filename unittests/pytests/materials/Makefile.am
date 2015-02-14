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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testmaterials.py
dist_check_SCRIPTS = testmaterials.py

noinst_PYTHON = \
	TestHomogeneous.py \
	TestMaterial.py \
	TestElasticMaterial.py \
	TestElasticPlaneStrain.py \
	TestElasticPlaneStress.py \
	TestElasticIsotropic3D.py \
	TestMaxwellIsotropic3D.py \
	TestMaxwellPlaneStrain.py \
	TestGenMaxwellIsotropic3D.py \
	TestGenMaxwellPlaneStrain.py \
	TestGenMaxwellQpQsIsotropic3D.py \
	TestPowerLaw3D.py \
	TestPowerLawPlaneStrain.py \
	TestDruckerPrager3D.py \
	TestDruckerPragerPlaneStrain.py


# End of file 
