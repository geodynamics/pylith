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
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

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
	shearnoslip_tri.cfg


noinst_TMP =


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
