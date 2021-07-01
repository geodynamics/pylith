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
	tri_onecell.mesh \
	tri_fourcells.mesh \
	tri.mesh \
	quad_onecell.mesh \
	quad_fourcells.mesh \
	quad_small.mesh \
	quad_aligned.mesh \
	tet_small.mesh \
	hex.mesh

noinst_DATA =

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/libtests/materials/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file
