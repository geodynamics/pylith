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
	elasticstrain1d.spatialdb \
	elasticplanestrain.spatialdb \
	elasticisotropic3d.spatialdb

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/feassemble/data
export-data: $(dist_noinst_DATA)
	for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(dist_noinst_DATA)


# End of file 
