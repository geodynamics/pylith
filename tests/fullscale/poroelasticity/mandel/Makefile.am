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
	generate_gmsh.py \
	meshes.py \
	TestMandel.py \
	mandel_soln.py \
	mandel_gendb.py

dist_noinst_DATA = \
	pylithapp.cfg \
	mesh_tri.msh \
	mesh_quad.msh \
	mandel_tri.cfg \
	mandel_quad.cfg

noinst_TMP = \
	mandel_disp.timedb


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
