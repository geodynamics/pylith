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
	TestThreeBlocks.py \
	threeblocks_soln.py \
	TestShearNoSlip.py \
	shearnoslip_soln.py \
	shearnoslip_gendb.py

dist_noinst_DATA = \
	mesh_tri.msh \
	mesh_quad.msh \
	geometry.jou \
	bc.jou \
	mesh_tri.jou \
	mesh_tri.exo \
	mesh_quad.jou \
	mesh_quad.exo \
	pylithapp.cfg \
	output_points.txt \
	threeblocks.cfg \
	threeblocks_quad.cfg \
	threeblocks_tri.cfg \
	threeblocks_cubit_quad.cfg \
	threeblocks_cubit_tri.cfg \
	threeblocks_ic.cfg \
	threeblocks_ic_quad.cfg \
	threeblocks_ic_tri.cfg \
	shearnoslip.cfg \
	shearnoslip_quad.cfg \
	shearnoslip_tri.cfg \
	zeroslipfn.timedb


noinst_TMP =


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
