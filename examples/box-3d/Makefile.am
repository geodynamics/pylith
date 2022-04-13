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
	bc.jou \
	mesh_hex.exo \
	mesh_hex.jou \
	mesh_tet.exo \
	mesh_tet.jou \
	pylithapp.cfg \
	step01_axialdisp.cfg \
	step02_sheardisp.cfg \
	step03_sheardisptract.cfg \
	step04_sheardispic.cfg \
	step05_sheardisptractrate.cfg \
	sheardisp_bc_xneg.spatialdb \
	sheardisp_bc_xpos.spatialdb \
	sheardisp_bc_yneg.spatialdb \
	sheardisp_bc_ypos.spatialdb \
	sheardisp_ic.spatialdb \
	sheardisprate_bc_xneg.spatialdb \
	sheardisprate_bc_yneg.spatialdb \
	viz/plot_dispwarp.py


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/examples/3d/box
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
