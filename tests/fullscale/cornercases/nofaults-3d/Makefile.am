# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	TestAxialDispOneCell.py \
	TestAxialDispConstrained.py \
	axialdisp_soln.py \
	axialdisp_gendb.py

dist_noinst_DATA = \
	fivecells_tet.mesh \
	onecell_hex.mesh \
	pylithapp.cfg \
	axialdisp.cfg \
	axialdisp_hex.cfg \
	dofconstrained.cfg \
	dofconstrained_tet.cfg \
	dofconstrained_hex.cfg

noinst_TMP = \
	axialdisp_bc.spatialdb



export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
