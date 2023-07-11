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
	TestTerzaghi.py \
	terzaghi_soln.py \
	terzaghi_gendb.py \
	TestTerzaghiCompaction.py \
	terzaghi_compaction_soln.py \
	terzaghi_compaction_gendb.py	

dist_noinst_DATA = \
	geometry.jou \
	bc.jou \
	mesh_tri.jou \
	mesh_tri.exo \
	mesh_quad.jou \
	mesh_quad.exo \
	terzaghi.cfg \
	terzaghi_tri.cfg \
	terzaghi_quad.cfg \
	terzaghi_compaction.cfg \
	terzaghi_compaction_tri.cfg \
	terzaghi_compaction_quad.cfg	

noinst_TMP = \
	terzaghi_bc.spatialdb \
	terzaghi_ic.spatialdb \
	terzaghi_compaction_bc.spatialdb \
	terzaghi_compaction_ic.spatialdb	



export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
