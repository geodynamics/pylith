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
