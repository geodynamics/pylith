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
	TestBodyForce.py \
	TestGravity.py \
	TestGravityBodyForce.py \
	bodyforce_soln.py \
	gravity_soln.py \
	gravity_bodyforce_soln.py

dist_noinst_DATA = \
	mesh_tri.msh \
	mesh_quad.msh \
	pylithapp.cfg \
	bodyforce.cfg \
	bodyforce_quad.cfg \
	bodyforce_tri.cfg \
	gravity.cfg \
	gravity_quad.cfg \
	gravity_tri.cfg \
	gravity_bodyforce.cfg \
	gravity_bodyforce_quad.cfg \
	gravity_bodyforce_tri.cfg


noinst_TMP =


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
