# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

dist_noinst_DATA = \
	cube2_ascii.gmv \
	cube2_ascii.pset \
	cube2_binary.gmv \
	cube2_binary.pset \
	cube2_binary_32on64.gmv \
	cube2_binary_32on64.pset \
	cube2_binary_64.gmv \
	cube2_binary_64.pset \
	mesh3DIndexOne.txt \
	twotri3.exo \
	twoquad4.exo \
	twotet4.exo \
	twohex8.exo \
	mesh2D_comments.txt \
	line2.mesh \
	line2_t10.vtk \
	line2_vertex_t10.vtk \
	line2_cell_t10.vtk \
	tri3.mesh \
	tri3_t10.vtk \
	tri3_vertex_t10.vtk \
	tri3_cell_t10.vtk \
	quad4.mesh \
	quad4_t10.vtk \
	quad4_vertex_t10.vtk \
	quad4_cell_t10.vtk \
	tet4.mesh \
	tet4_t10.vtk \
	tet4_vertex_t10.vtk \
	tet4_cell_t10.vtk \
	hex8.mesh \
	hex8_t10.vtk \
	hex8_vertex_t10.vtk \
	hex8_cell_t10.vtk \
	line2_mat_t10.vtk \
	line2_mat_vertex_t10.vtk \
	line2_mat_cell_t10.vtk \
	tri3_mat_t10.vtk \
	tri3_mat_vertex_t10.vtk \
	tri3_mat_cell_t10.vtk \
	quad4_mat_t10.vtk \
	quad4_mat_vertex_t10.vtk \
	quad4_mat_cell_t10.vtk \
	tet4_mat_t10.vtk \
	tet4_mat_vertex_t10.vtk \
	tet4_mat_cell_t10.vtk \
	hex8_mat_t10.vtk \
	hex8_mat_vertex_t10.vtk \
	hex8_mat_cell_t10.vtk \
	tri3_fault_t10.vtk \
	tri3_fault_vertex_t10.vtk \
	tri3_fault_cell_t10.vtk \
	quad4_fault_t10.vtk \
	quad4_fault_vertex_t10.vtk \
	quad4_fault_cell_t10.vtk \
	tet4_fault_t10.vtk \
	tet4_fault_vertex_t10.vtk \
	tet4_fault_cell_t10.vtk \
	hex8_fault_t10.vtk \
	hex8_fault_vertex_t10.vtk \
	hex8_fault_cell_t10.vtk \
	tri3_bc_t10.vtk \
	tri3_bc_vertex_t10.vtk \
	tri3_bc_cell_t10.vtk \
	quad4_bc_t10.vtk \
	quad4_bc_vertex_t10.vtk \
	quad4_bc_cell_t10.vtk \
	tet4_bc_t10.vtk \
	tet4_bc_vertex_t10.vtk \
	tet4_bc_cell_t10.vtk \
	hex8_bc_t10.vtk \
	hex8_bc_vertex_t10.vtk \
	hex8_bc_cell_t10.vtk \
	output_t12.vtk \
	output_vertex_t12.vtk \
	output_vertex_filter_t12.vtk \
	output_cell_t12.vtk \
	output_cell_filter_t12.vtk

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/meshio/data
export-data: $(dist_noinst_DATA)
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file 
