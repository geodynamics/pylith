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
	mesh_tri.msh \
	mesh_quad.msh \
	geometry.jou \
	bc.jou \
	mesh_tri.jou \
	mesh_tri.exo \
	mesh_quad.jou \
	mesh_quad.exo \
	pylithapp.cfg \
	axialdisp.cfg \
	axialdisp_tri.cfg \
	axialdisp_quad.cfg \
	axialdisp_cubit_tri.cfg \
	axialdisp_cubit_quad.cfg \
	sheartraction.cfg \
	sheartraction_tri.cfg \
	sheartraction_quad.cfg \
	sheartraction_rate.cfg \
	sheartraction_rate_tri.cfg \
	sheartraction_rate_quad.cfg \
	gravity.cfg \
	gravity_tri.cfg \
	gravity_quad.cfg \
	gravity_refstate.cfg \
	gravity_refstate_tri.cfg \
	gravity_refstate_quad.cfg \
	output_points.txt

noinst_TMP = \
	axialdisp_bc.spatialdb \
	axialdisp_ic.spatialdb \
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
