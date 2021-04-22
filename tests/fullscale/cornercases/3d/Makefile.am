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
	TestAxialDispRotated.py \
	axialdisp_soln.py \
	axialdisp_cylinder_gendb.py \
	axialdisp_pgram_gendb.py \
	axialdisp_rotated_soln.py \
	axialdisp_rotated_gendb.py

dist_noinst_DATA = \
	fivecells_tet.mesh \
	onecell_hex.mesh \
	fivecells_rotated_tet.mesh \
	onecell_rotated_hex.mesh \
	createbc_cylinder.jou \
	createbc_pgram.jou \
	geometry_cylinder.jou \
	geometry_pgram.jou \
	mesh_cylinder_hex.jou \
	mesh_cylinder_hex.exo \
	mesh_cylinder_tet.jou \
	mesh_cylinder_tet.exo \
	mesh_pgram_hex.jou \
	mesh_pgram_hex.exo \
	mesh_pgram_tet.jou \
	mesh_pgram_tet.exo \
	pylithapp.cfg \
	axialdisp.cfg \
	axialdisp_hex.cfg \
	axialdisp_cylinder.cfg \
	axialdisp_cylinder_hex.cfg \
	axialdisp_cylinder_tet.cfg \
	axialdisp_pgram.cfg \
	axialdisp_pgram_hex.cfg \
	axialdisp_pgram_tet.cfg \
	axialdisp_rotated.cfg \
	axialdisp_rotated_hex.cfg \
	axialdisp_rotated_tet.cfg \
	dofconstrained.cfg \
	dofconstrained_tet.cfg \
	dofconstrained_hex.cfg

noinst_TMP = \
	axialdisp_bc.spatialdb \
	axialdisp_cylinder_bc.spatialdb \
	axialdisp_pgram_bc.spatialdb \
	axialdisp_rotated_bc.spatialdb



# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/fullscale/cornercases/3d
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
