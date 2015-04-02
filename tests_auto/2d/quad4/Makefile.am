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

TESTS = testpylith.py

dist_check_SCRIPTS = testpylith.py

dist_noinst_PYTHON = \
	TestQuad4.py \
	TestAxialDisp.py \
	axialdisp_soln.py \
	axialdisp_gendb.py \
	TestShearDisp.py \
	sheardisp_soln.py \
	sheardisp_gendb.py \
	TestDislocation.py \
	TestDislocationTwoFaults.py \
	dislocation_soln.py \
	dislocation2_soln.py \
	TestLgDeformRigidBody.py \
	rigidbody_soln.py \
	rigidbody_gendb.py \
	TestLgDeformTraction.py \
	lgdeformtraction_soln.py \
	TestFrictionCompression.py \
	TestFrictionOpening.py \
	TestFrictionShearSliding.py \
	TestFrictionShearStick.py \
	TestSlipWeakeningCompression.py \
	TestSlipWeakeningOpening.py \
	TestSlipWeakeningShearSliding.py \
	TestSlipWeakeningShearStick.py \
	friction_compression_soln.py \
	friction_shear_stick_soln.py \
	friction_shear_sliding_soln.py \
	friction_opening_soln.py \
	slipweakening_compression_soln.py \
	slipweakening_shear_stick_soln.py \
	slipweakening_shear_sliding_soln.py \
	slipweakening_opening_soln.py

dist_noinst_DATA = \
	geometry.jou \
	mesh.jou \
	mesh.exo \
	matprops.spatialdb \
	axialdisp.cfg \
	sheardisp.cfg \
	dislocation.cfg \
	dislocation_np2.cfg \
	dislocation2.cfg \
	lgdeformrigidbody.cfg \
	lgdeformtraction.cfg \
	friction_compression.cfg \
	friction_shear_stick.cfg \
	friction_shear_sliding.cfg \
	friction_opening.cfg \
	slipweakening_compression.cfg \
	slipweakening_shear_stick.cfg \
	slipweakening_shear_sliding.cfg \
	slipweakening_opening.cfg


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests_auto/2d/quad4
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.h5 *.xmf *.pyc


# End of file 
