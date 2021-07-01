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
#

dist_noinst_DATA = \
	README.md \
	pylithapp.cfg \
	step01_coseismic.cfg \
	step02_interseismic.cfg \
	step03_eqcycle.cfg \
	step04_afterslip.cfg \
	step05_eqcycleslipweakening.cfg \
	step06_eqcycleratestate.cfg \
	geometry.jou \
	createbc.jou \
	mesh_tri.jou \
	mesh_tri.exo \
	fault_coseismic.spatialdb \
	fault_slabtop_creep.spatialdb \
	fault_slabtop_ratestate.spatialdb \
	fault_slabtop_slipweakening.spatialdb \
	fault_slabtop_tractions.spatialdb \
	mat_concrust.spatialdb \
	mat_conmantle.spatialdb \
	mat_oceancrust.spatialdb \
	mat_oceanmantle.spatialdb \
	afterslip_tractions.spatialdb \
	afterslip_tractions.py


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/examples/2d/subduction
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
