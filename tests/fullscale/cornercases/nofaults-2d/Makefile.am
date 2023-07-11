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



export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
