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
	injection_2d.txt \
	extraction_2d.txt \
	pylithapp.cfg \
	step00_steadystate.cfg \
	step02_5spot_injection_BC.cfg \
	step03_5spot_injection_BC_enhanced.cfg \
	mesh/geometry.jou \
	mesh/bc.jou \
	mesh/mesh_quad.jou \
	mesh/mesh_quad.exo \
	mesh/mesh_enhanced_quad.jou \
	mesh/mesh_enhanced_quad.exo \
    spatialdb/spe_perm.dat \
    spatialdb/spe_phi.dat \
    spatialdb/SPE10_grid_pylith.py \
    spatialdb/SPE10_pylith.py \
    spatialdb/SPE10_parameters_grid.spatialdb \
    spatialdb/SPE10_parameters.spatialdb


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/examples/SPE10-2d
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA)  $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.pyc output/*


# End of file 
