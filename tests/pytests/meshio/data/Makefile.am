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
	mesh2Din3D.txt \
	cube2_ascii.gmv \
	cube2_ascii.pset \
	cube2.txt \
	twohex8.exo \
	twohex8.txt \
	point.txt \
	points.txt \
	tri3.h5 \
	tri3.xmf \
	tri3_vertex.h5 \
	tri3_vertex.xmf \
	tri3_cell.h5 \
	tri3_cell.xmf \
	tri3_points.h5 \
	tri3_points.xmf \
	tri3_points_vertex.h5 \
	tri3_points_vertex.xmf \
	tri3_surf.h5 \
	tri3_surf.xmf \
	tri3_surf_vertex.h5 \
	tri3_surf_vertex.xmf \
	tri3_surf_cell.h5 \
	tri3_surf_cell.xmf \
	quad4.h5 \
	quad4.xmf \
	quad4_vertex.h5 \
	quad4_vertex.xmf \
	quad4_cell.h5 \
	quad4_cell.xmf \
	quad4_points.h5 \
	quad4_points.xmf \
	quad4_points_vertex.h5 \
	quad4_points_vertex.xmf \
	quad4_surf.h5 \
	quad4_surf.xmf \
	quad4_surf_vertex.h5 \
	quad4_surf_vertex.xmf \
	quad4_surf_cell.h5 \
	quad4_surf_cell.xmf \
	tet4.h5 \
	tet4.xmf \
	tet4_vertex.h5 \
	tet4_vertex.xmf \
	tet4_cell.h5 \
	tet4_cell.xmf \
	tet4_points.h5 \
	tet4_points.xmf \
	tet4_points_vertex.h5 \
	tet4_points_vertex.xmf \
	tet4_surf.h5 \
	tet4_surf.xmf \
	tet4_surf_vertex.h5 \
	tet4_surf_vertex.xmf \
	tet4_surf_cell.h5 \
	tet4_surf_cell.xmf \
	hex8.h5 \
	hex8.xmf \
	hex8_vertex.h5 \
	hex8_vertex.xmf \
	hex8_cell.h5 \
	hex8_cell.xmf \
	hex8_points.h5 \
	hex8_points.xmf \
	hex8_points_vertex.h5 \
	hex8_points_vertex.xmf \
	hex8_surf.h5 \
	hex8_surf.xmf \
	hex8_surf_vertex.h5 \
	hex8_surf_vertex.xmf \
	hex8_surf_cell.h5 \
	hex8_surf_cell.xmf

noinst_TMP = \
	cube2_test.txt \
	mesh2Din3D_test.txt \
	twohex8_test.txt


# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/pytests/meshio/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file
