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
	TestMandel.py \
	mandel_soln.py \
	mandel_gendb.py \
	TestMandelCompaction.py \
	mandel_compaction_soln.py \
	mandel_compaction_gendb.py	

dist_noinst_DATA = \
	geometry.jou \
	bc.jou \
	mesh_tri.jou \
	mesh_tri.exo \
	mesh_quad.jou \
	mesh_quad.exo \
	mandel.cfg \
	mandel_tri.cfg \
	mandel_quad.cfg \
	mandel_compaction.cfg \
	mandel_compaction_tri.cfg \
	mandel_compaction_quad.cfg 	

noinst_TMP = \
	mandel_bc.spatialdb \
	mandel_ic.spatialdb \
	mandel_compaction_bc.spatialdb \
	mandel_compaction_ic.spatialdb	



export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
