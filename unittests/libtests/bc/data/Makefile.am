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
	line2.mesh \
	line2_disp.spatialdb \
	line2_tractions.spatialdb \
	elasticstrain1d.spatialdb \
	tri3.mesh \
	tri3_disp.spatialdb \
	tri3_vel.spatialdb \
	tri3_disp2.spatialdb \
	tri3_vel2.spatialdb \
	tri3_tractions.spatialdb \
	tri3_force.spatialdb \
	elasticplanestrain.spatialdb \
	quad4.mesh \
	quad4_disp.spatialdb \
	quad4_tractions.spatialdb \
	tet4.mesh \
	tet4_disp.spatialdb \
	tet4_disp2.spatialdb \
	tet4_vel2.spatialdb \
	tet4_tractions.spatialdb \
	elasticisotropic3d.spatialdb \
	hex8.mesh \
	hex8_disp.spatialdb \
	hex8b.mesh \
	hex8b_disp.spatialdb


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/bc/data
export-data: $(dist_noinst_DATA)
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file 
