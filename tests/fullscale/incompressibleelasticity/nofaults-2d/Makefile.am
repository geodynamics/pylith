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
	TestGravityIncompressible.py \
	gravity_incompressible_soln.py

dist_noinst_DATA = \
	mesh_tri.msh \
	mesh_quad.msh \
	pylithapp.cfg \
	gravity_incompressible.cfg \
	gravity_incompressible_quad.cfg \
	gravity_incompressible_tri.cfg \
	gravity_incompressible_ic.cfg \
	gravity_incompressible_ic_quad.cfg \
	gravity_incompressible_ic_tri.cfg \
	initial_conditions.spatialdb \
	output_points.txt

noinst_TMP = 


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
