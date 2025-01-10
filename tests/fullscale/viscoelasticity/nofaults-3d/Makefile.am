# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	TestAxialTractionMaxwell.py \
	axialtraction_maxwell_soln.py \
	axialtraction_maxwell_gendb.py \
	TestAxialStrainGenMaxwell.py \
	axialstrain_genmaxwell_soln.py \
	axialstrain_genmaxwell_gendb.py \
	TestAxialStrainRateGenMaxwell.py \
	axialstrainrate_genmaxwell_soln.py \
	axialstrainrate_genmaxwell_gendb.py

dist_noinst_DATA = \
	geometry.jou \
	bc.jou \
	mesh_tet.jou \
	mesh_tet.exo \
	mesh_hex.jou \
	mesh_hex.exo \
	pylithapp.cfg \
	axialtraction_maxwell.cfg \
	axialtraction_maxwell_tet.cfg \
	axialtraction_maxwell_hex.cfg \
	mat_maxwell.spatialdb \
	axialstrain_genmaxwell.cfg \
	axialstrain_genmaxwell_tet.cfg \
	axialstrain_genmaxwell_hex.cfg \
	axialstrainrate_genmaxwell.cfg \
	axialstrainrate_genmaxwell_tet.cfg \
	axialstrainrate_genmaxwell_hex.cfg \
	mat_genmaxwell.spatialdb

noinst_TMP = \
	axialtraction_maxwell_bc.spatialdb \
	axialstrain_genmaxwell_bc.spatialdb \
	axialstrainrate_genmaxwell_bc.spatialdb



export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
