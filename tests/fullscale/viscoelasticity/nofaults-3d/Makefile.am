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



# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/fullscale/viscoelasticity/nofaults-3d
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
