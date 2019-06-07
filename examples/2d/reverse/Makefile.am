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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

include $(top_srcdir)/check.am

TESTS = run_examples.py

dist_check_SCRIPTS = test_pylith.py

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
	step01_gravity.cfg \
	step02_gravity_refstate.cfg \
	step03_gravity_incompressible.cfg \
	step04_neumann.cfg \
	step05_onefault.cfg \
	step06_twofault.cfg \
	elastic_incompressible_matfields.spatialdb \
	elastic_matfields.spatialdb \
	gravity_refstate_matfields.spatialdb \
	traction_trapezoid.spatialdb \
	viz/plot_dispwarp.py


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/examples/2d/box
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
