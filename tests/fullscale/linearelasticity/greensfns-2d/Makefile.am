# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
