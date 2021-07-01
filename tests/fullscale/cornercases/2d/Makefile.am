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
	TestAxialDispOneCell.py \
	TestAxialDispConstrained.py \
	axialdisp_soln.py \
	axialdisp_gendb.py

dist_noinst_DATA = \
	twocells_tri.mesh \
	onecell_quad.mesh \
	pylithapp.cfg \
	axialdisp.cfg \
	axialdisp_quad.cfg \
	dofconstrained.cfg \
	dofconstrained_tri.cfg \
	dofconstrained_quad.cfg

noinst_TMP = \
	axialdisp_bc.spatialdb



# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/fullscale/cornercases/2d
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
