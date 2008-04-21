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
	bulkprops_2d.spatialdb \
	tri3.mesh \
	tri3_finalslip.spatialdb \
	tri3_sliptime.spatialdb \
	tri3_peakrate.spatialdb

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/pytests/faults/data
export-data: $(dist_noinst_DATA)
	for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(dist_noinst_DATA)


# End of file 
