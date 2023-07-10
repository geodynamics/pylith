# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	generate_gmsh.py \
	meshes.py \
	TestLeftLateral.py \
	TestOpening.py \
	TestSlipThreshold.py \
	faultimpulses_soln.py

dist_noinst_DATA = \
	mesh_tri.msh \
	mesh_quad.msh \
	pylithapp.cfg \
	leftlateral_b0.cfg \
	leftlateral_b0_quad.cfg \
	leftlateral_b0_tri.cfg \
	leftlateral_b1.cfg \
	leftlateral_b1_quad.cfg \
	leftlateral_b1_tri.cfg \
	leftlateral_b2.cfg \
	leftlateral_b2_quad.cfg \
	leftlateral_b2_tri.cfg \
	opening.cfg \
	opening_quad.cfg \
	opening_tri.cfg \
	slipthreshold.cfg \
	slipthreshold_quad.cfg \
	slipthreshold_tri.cfg \
	slip_ypos.spatialdb


noinst_TMP =


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
