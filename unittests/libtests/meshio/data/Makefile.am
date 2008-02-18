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

noinst_DATA = \
	cube2_ascii.gmv \
	cube2_ascii.pset \
	cube2_binary.gmv \
	cube2_binary.pset \
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
	hex8_cell_t10.vtk

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/meshio/data
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA)


# End of file 
