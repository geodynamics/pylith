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
	line2.mesh \
	line2.spatialdb \
	tri3.mesh \
	tri3.spatialdb \
	tri3_b.spatialdb \
	quad4.mesh \
	quad4.spatialdb \
	tet4.mesh \
	tet4.spatialdb \
	hex8.mesh \
	hex8.spatialdb

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/bc/data
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA)


# End of file 
