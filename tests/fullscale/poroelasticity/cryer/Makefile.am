# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	TestCryer.py \
	cryer_soln.py 

dist_noinst_DATA = \
	geometry.jou \
	bc.jou \
	mesh_tet.jou \
	mesh_tet.exo \
	mesh_hex.jou \
	mesh_hex.exo \
	cryer.cfg \
	cryer_tet.cfg \
	cryer_hex.cfg

noinst_TMP = \
	cryer_bc.spatialdb \
	cryer_ic.spatialdb

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
