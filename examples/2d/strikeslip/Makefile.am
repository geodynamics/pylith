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
	mesh_tri.jou \
	mesh_quad.jou \
	gradient.jou \
	createbc.jou \
	pylithapp.cfg \
	step01_slip.cfg \
	step02_slip_velbc.cfg \
	step03_multislip_velbc.cfg \
	disprate_bc_xneg.spatialdb \
	disprate_bc_xpos.spatialdb \
	elastic_xneg_matfields.spatialdb \
	elastic_xpos_matfields.spatialdb \
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