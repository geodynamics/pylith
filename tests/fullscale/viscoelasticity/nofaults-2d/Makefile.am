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
	meshes.py \
	TestAxialTractionMaxwell.py \
	TestAxialStrainGenMaxwell.py \
	TestAxialStrainRateGenMaxwell.py \
	axialtraction_maxwell_soln.py \
	axialtraction_maxwell_gendb.py \
	axialstrain_genmaxwell_soln.py \
	axialstrain_genmaxwell_gendb.py \
	axialstrainrate_genmaxwell_soln.py \
	axialstrainrate_genmaxwell_gendb.py

dist_noinst_DATA = \
	geometry.jou \
	bc.jou \
	mesh_tri.jou \
	mesh_tri.exo \
	mesh_quad.jou \
	mesh_quad.exo \
	pylithapp.cfg \
	axialtraction_maxwell.cfg \
	axialtraction_maxwell_tri.cfg \
	axialtraction_maxwell_quad.cfg \
	axialstrain_genmaxwell.cfg \
	axialstrain_genmaxwell_tri.cfg \
	axialstrain_genmaxwell_quad.cfg \
	axialstrainrate_genmaxwell.cfg \
	axialstrainrate_genmaxwell_tri.cfg \
	axialstrainrate_genmaxwell_quad.cfg \
	mat_maxwell.spatialdb \
	mat_genmaxwell.spatialdb

noinst_TMP = \
	axialtraction_maxwell_disp.spatialdb \
	axialstrain_genmaxwell_disp.spatialdb \
	axialstrain_genmaxwell_bc.spatialdb \
	axialstrainrate_genmaxwell_disp.spatialdb \
	axialstrainrate_genmaxwell_bc.spatialdb



export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
