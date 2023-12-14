# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

TESTS = test_eqinfo.py

dist_check_SCRIPTS = test_eqinfo.py

dist_noinst_PYTHON = \
	TestEqInfo.py \
	TestEqInfoLine.py \
	TestEqInfoTri.py \
	TestEqInfoQuad.py

dist_noinst_DATA = \
	mat_elastic.spatialdb \
	line.cfg \
	line_one.h5 \
	line_two.h5 \
	tri.cfg \
	tri_one.h5 \
	tri_two.h5 \
	quad.cfg \
	quad_one.h5 \
	quad_two.h5

noinst_TMP = \
	stats_line.py \
	stats_tri.py \
	stats_quad.py


TESTS_ENVIRONMENT = PYTHONPATH='$(PYTHONPATH):$(abs_builddir)'

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
