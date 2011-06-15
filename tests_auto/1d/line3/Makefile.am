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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

TESTS = testpylith.py

check_SCRIPTS = testpylith.py

dist_noinst_PYTHON = \
	TestLine3.py \
	TestAxial.py \
	TestDislocation.py \
	testpylith.py

dist_noinst_DATA = \
	bar.mesh \
	matprops.spatialdb \
	axialextension.cfg \
	axialextension_disp.spatialdb \
	dislocation.cfg \
	finalslip.spatialdb \
	sliptime.spatialdb

noinst_TMP = \
	axialextension_t0000000.vtk \
	axialextension-elastic_info.vtk \
	axialextension-elastic_t0000000.vtk \
	dislocation-fault_info.vtk \
	dislocation-fault_t0000000.vtk \
	dislocation-elastic_info.vtk \
	dislocation-elastic_t0000000.vtk \
	dislocation_t0000000.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests_auto/1d/line3
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data
CLEANFILES = *.pyc


# End of file 
