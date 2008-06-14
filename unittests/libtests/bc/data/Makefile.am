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
	line2.spatialdb \
	line2-tractions.spatialdb \
	elasticstrain1d.spatialdb \
	tri3.mesh \
	tri3.spatialdb \
	tri3_rate.spatialdb \
	tri3_b.spatialdb \
	tri3_b_rate.spatialdb \
	tri3-tractions.spatialdb \
	elasticplanestrain.spatialdb \
	quad4.mesh \
	quad4.spatialdb \
	quad4-tractions.spatialdb \
	tet4.mesh \
	tet4.spatialdb \
	tet4_a.spatialdb \
	tet4_a_rate.spatialdb \
	tet4-tractions.spatialdb \
	elasticisotropic3d.spatialdb \
	hex8.mesh \
	hex8.spatialdb \
	hex8b.mesh \
	hex8b.spatialdb

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/bc/data
export-data: $(dist_noinst_DATA)
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file 
