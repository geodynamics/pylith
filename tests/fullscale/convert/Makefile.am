# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

TESTS = test_convert.py

dist_check_SCRIPTS = test_convert.py

dist_noinst_PYTHON = \
	TestConvert.py

dist_noinst_DATA = \
	convert_quad.cfg \
	quad_small.mesh

noinst_TMP = \
	quad_small.h5


TESTS_ENVIRONMENT = PYTHONPATH='$(PYTHONPATH):$(abs_builddir)'

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
