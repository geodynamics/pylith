# -*- Makefile -*-
#
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
#

dist_noinst_DATA = \
	bulkprops_2d.spatialdb \
	tri3.mesh \
	tri3_finalslip.spatialdb \
	tri3_sliptime.spatialdb \
	tri3_peakrate.spatialdb \
	tri3_initialtractions.spatialdb \
	tri3_changetractions.spatialdb \
	tri3_staticfriction.spatialdb \
	tri3_impulses.spatialdb \
	slipfn.timedb

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/pytests/faults/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file 
