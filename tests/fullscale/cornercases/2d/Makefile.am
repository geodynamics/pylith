# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2019 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

include $(top_srcdir)/check.am

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	TestAxialDispOneCell.py \
	TestAxialDispConstrained.py \
	TestAxialDispCylinder.py \
	TestAxialDispPgram.py \
	axialdisp_soln.py \
	axialdisp_gendb.py \
	axialdisp_cylinder_gendb.py \
	axialdisp_pgram_gendb.py

dist_noinst_DATA = \
	twocells_tri.mesh \
	onecell_quad.mesh \
	createbc_cylinder.jou \
	createbc_pgram.jou \
	geometry_cylinder.jou \
	geometry_pgram.jou \
	mesh_cylinder_quad.jou \
	mesh_cylinder_tri.jou \
	mesh_cylinder_quad.exo \
	mesh_cylinder_tri.exo \
	mesh_pgram_quad.jou \
	mesh_pgram_tri.jou \
	mesh_4cells_pgram_quad.exo \
	mesh_8cells_pgram_tri.exo \
	pylithapp.cfg \
	axialdisp.cfg \
	axialdisp_quad.cfg \
	axialdisp_cylinder.cfg \
	axialdisp_cylinder_quad.cfg \
	axialdisp_cylinder_tri.cfg \
	axialdisp_pgram.cfg \
	axialdisp_pgram_4cells_quad.cfg \
	axialdisp_pgram_8cells_tri.cfg \
	dofconstrained.cfg \
	dofconstrained_tri.cfg \
	dofconstrained_quad.cfg

noinst_TMP = \
	axialdisp_bc.spatialdb \
	axialdisp_cylinder_bc.spatialdb \
	axialdisp_pgram_bc.spatialdb



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
