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
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

dist_noinst_DATA = \
	README.md \
	geometry.jou \
	gradient.jou \
	createbc.jou \
	mesh_quad.jou \
	mesh_quad.exo \
	mesh_tri.exo \
	mesh_tri.jou \
	pylithapp.cfg \
	solver_onefault.cfg \
	solver_twofaults.cfg \
	step01_gravity.cfg \
	step02_gravity_refstate.cfg \
	step03_gravity_incompressible.cfg \
	step04_surfload.cfg \
	step05_onefault.cfg \
	step06_twofaults_elastic.cfg \
	step07_twofaults_maxwell.cfg \
	step08_twofaults_powerlaw.cfg \
	mat_elastic_incompressible.spatialdb \
	mat_elastic.spatialdb \
	mat_gravity_refstate.spatialdb \
	mat_maxwell.spatialdb \
	mat_powerlaw.spatialdb \
	traction_surfload.spatialdb \
	viz/plot_dispwarp.py \
	viz/plot_mesh.py


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/examples/2d/reverse
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
