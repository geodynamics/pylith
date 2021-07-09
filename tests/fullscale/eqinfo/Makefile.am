# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

include $(top_srcdir)/check.am

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

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/fullscale/eqinfo
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data
clean-local: clean-data
CLEANFILES = *.pyc

# End of file
