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
	line2_finalslip.spatialdb \
	line2_sliptime.spatialdb \
	line2_peakrate.spatialdb \
	tri3.mesh \
	tri3_finalslip.spatialdb \
	tri3_sliptime.spatialdb \
	tri3_peakrate.spatialdb \
	quad4.mesh \
	tet4.mesh \
	tet4_finalslip.spatialdb \
	tet4_sliptime.spatialdb \
	tet4_peakrate.spatialdb \
	hex8.mesh

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/faults/data
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
