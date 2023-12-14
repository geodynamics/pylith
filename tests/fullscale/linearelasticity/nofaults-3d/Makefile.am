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
	meshes.py \
	TestAxialDisp.py \
	axialdisp_soln.py \
	axialdisp_gendb.py \
	TestShearTraction.py \
	sheartraction_soln.py \
	sheartraction_gendb.py \
	TestShearTractionRate.py \
	sheartraction_rate_soln.py \
	sheartraction_rate_gendb.py \
	TestGravity.py \
	gravity_soln.py \
	TestGravityRefState.py \
	gravity_refstate_soln.py \
	gravity_refstate_gendb.py

dist_noinst_DATA = \
	geometry.jou \
	bc.jou \
	mesh_tet.jou \
	mesh_tet.exo \
	mesh_hex.jou \
	mesh_hex.exo \
	output_points.txt \
	pylithapp.cfg \
	axialdisp.cfg \
	axialdisp_tet.cfg \
	axialdisp_hex.cfg \
	sheartraction.cfg \
	sheartraction_tet.cfg \
	sheartraction_hex.cfg \
 	sheartraction_rate.cfg \
	sheartraction_rate_tet.cfg \
	sheartraction_rate_hex.cfg \
	gravity.cfg \
	gravity_tet.cfg \
	gravity_hex.cfg \
	gravity_refstate.cfg \
	gravity_refstate_tet.cfg \
	gravity_refstate_hex.cfg

noinst_TMP = \
	axialdisp_bc.spatialdb \
	axialdisp_ic.spatialdb
	sheartraction_disp.spatialdb \
	sheartraction_rate_disp.spatialdb \
	gravity_refstate_matfields.spatialdb



export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
