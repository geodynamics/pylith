# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

include $(top_srcdir)/check.am

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
	gravity_refstate_gendb.py \
	TestGravityIncompressible.py \
	gravity_incompressible_soln.py

dist_noinst_DATA = \
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
	gravity_incompressible.cfg \
	gravity_incompressible_quad.cfg \
	gravity_incompressible_tri.cfg \
	output_points.txt

noinst_TMP = \
	axialdisp_bc.spatialdb \
	axialdisp_ic.spatialdb \
	sheartraction_disp.spatialdb \
	sheartraction_rate_disp.spatialdb \
	gravity_refstate_matfields.spatialdb



# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/fullscale/linearelasticity/nofaults-2d
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.pyc output/*


# End of file
