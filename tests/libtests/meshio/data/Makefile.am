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
	mesh2D_comments.txt \
	mesh3D_index1.txt \
	cube2_ascii.gmv \
	cube2_ascii.pset \
	cube2_binary.gmv \
	cube2_binary.pset \
	cube2_binary_32on64.gmv \
	cube2_binary_32on64.pset \
	cube2_binary_64.gmv \
	cube2_binary_64.pset \
	twotri3_12.2.exo \
	twotri3_13.0.exo \
	twoquad4_12.2.exo \
	twoquad4_13.0.exo \
	twotet4_12.2.exo \
	twotet4_13.0.exo \
	twohex8_12.2.exo \
	twohex8_13.0.exo \
	mesh_tri3.exo \
	mesh_quad4.exo \
	mesh_tet4.exo \
	mesh_hex8.exo \
	tri3.mesh \
	tri3_vertex_t10.vtk \
	tri3_cell_t10.vtk \
	quad4.mesh \
	quad4_vertex_t10.vtk \
	quad4_cell_t10.vtk \
	tet4.mesh \
	tet4_vertex_t10.vtk \
	tet4_cell_t10.vtk \
	hex8.mesh \
	hex8_vertex_t10.vtk \
	hex8_cell_t10.vtk \
	tri3_mat_vertex_t10.vtk \
	tri3_mat_cell_t10.vtk \
	quad4_mat_vertex_t10.vtk \
	quad4_mat_cell_t10.vtk \
	tet4_mat_vertex_t10.vtk \
	tet4_mat_cell_t10.vtk \
	hex8_mat_vertex_t10.vtk \
	hex8_mat_cell_t10.vtk \
	tri3_surf_vertex_t10.vtk \
	tri3_surf_cell_t10.vtk \
	quad4_surf_vertex_t10.vtk \
	quad4_surf_cell_t10.vtk \
	tet4_surf_vertex_t10.vtk \
	tet4_surf_cell_t10.vtk \
	hex8_surf_vertex_t10.vtk \
	hex8_surf_cell_t10.vtk \
	tri3_bc_vertex_t10.vtk \
	tri3_bc_cell_t10.vtk \
	quad4_bc_vertex_t10.vtk \
	quad4_bc_cell_t10.vtk \
	tet4_bc_vertex_t10.vtk \
	tet4_bc_cell_t10.vtk \
	hex8_bc_vertex_t10.vtk \
	hex8_bc_cell_t10.vtk \
	tri3_fault_vertex_t10.vtk \
	tri3_fault_cell_t10.vtk \
	quad4_fault_vertex_t10.vtk \
	quad4_fault_cell_t10.vtk \
	tet4_fault_vertex_t10.vtk \
	tet4_fault_cell_t10.vtk \
	hex8_fault_vertex_t10.vtk \
	hex8_fault_cell_t10.vtk \
	output_t12.vtk \
	output_vertex_t12.vtk \
	output_vertex_filter_t12.vtk \
	output_cell_t12.vtk \
	output_cell_filter_t12.vtk \
	tri3.h5 \
	tri3_cell.h5 \
	tri3_vertex.h5 \
	quad4.h5 \
	quad4_cell.h5 \
	quad4_vertex.h5 \
	tet4.h5 \
	tet4_cell.h5 \
	tet4_vertex.h5 \
	hex8.h5 \
	hex8_cell.h5 \
	hex8_vertex.h5 \
	tri3_mat.h5 \
	tri3_mat_vertex.h5 \
	tri3_mat_cell.h5 \
	quad4_mat.h5 \
	quad4_mat_vertex.h5 \
	quad4_mat_cell.h5 \
	tet4_mat.h5 \
	tet4_mat_vertex.h5 \
	tet4_mat_cell.h5 \
	hex8_mat.h5 \
	hex8_mat_vertex.h5 \
	hex8_mat_cell.h5 \
	tri3_surf.h5 \
	tri3_surf_vertex.h5 \
	tri3_surf_cell.h5 \
	quad4_surf.h5 \
	quad4_surf_vertex.h5 \
	quad4_surf_cell.h5 \
	tet4_surf.h5 \
	tet4_surf_vertex.h5 \
	tet4_surf_cell.h5 \
	hex8_surf.h5 \
	hex8_surf_vertex.h5 \
	hex8_surf_cell.h5 \
	tri3_bc.h5 \
	tri3_bc_vertex.h5 \
	tri3_bc_cell.h5 \
	quad4_bc.h5 \
	quad4_bc_vertex.h5 \
	quad4_bc_cell.h5 \
	tet4_bc.h5 \
	tet4_bc_vertex.h5 \
	tet4_bc_cell.h5 \
	hex8_bc.h5 \
	hex8_bc_vertex.h5 \
	hex8_bc_cell.h5 \
	tri3_fault.h5 \
	tri3_fault_vertex.h5 \
	tri3_fault_cell.h5 \
	quad4_fault.h5 \
	quad4_fault_vertex.h5 \
	quad4_fault_cell.h5 \
	tet4_fault.h5 \
	tet4_fault_vertex.h5 \
	tet4_fault_cell.h5 \
	hex8_fault.h5 \
	hex8_fault_vertex.h5 \
	hex8_fault_cell.h5 \
	tri3_points_vertex_t10.vtk \
	quad4_points_vertex_t10.vtk \
	tet4_points_vertex_t10.vtk \
	hex8_points_vertex_t10.vtk \
	tri3_points.h5 \
	tri3_points_vertex.h5 \
	quad4_points.h5 \
	quad4_points_vertex.h5 \
	tet4_points.h5 \
	tet4_points_vertex.h5 \
	hex8_points.h5 \
	hex8_points_vertex.h5


noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/libtests/meshio/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file 
