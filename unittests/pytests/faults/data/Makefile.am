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
	tri3.mesh \
	tri3_finalslip.spatialdb \
	tri3_sliptime.spatialdb \
	tri3_peakrate.spatialdb

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/pytests/faults/data
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA)


# End of file 
