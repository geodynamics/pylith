# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

if ENABLE_HDF5
TESTS = test_eqinfo.py
endif

dist_check_SCRIPTS = test_eqinfo.py

dist_noinst_PYTHON = \
	TestEqInfo.py \
	TestEqInfoLine.py \
	TestEqInfoTri3.py \
	TestEqInfoQuad4.py

dist_noinst_DATA = \
	mat_elastic.spatialdb \
	line.cfg \
	line_one.h5 \
	line_two.h5 \
	tri3.cfg \
	tri3_one.h5 \
	tri3_two.h5 \
	quad4.cfg \
	quad4_one.h5 \
	quad4_two.h5

noinst_TMP = \
	stats_line.py \
	stats_tri3.py \
	stats_quad4.py


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests_auto/eqinfo
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data
clean-local: clean-data
CLEANFILES = *.pyc

# End of file 
