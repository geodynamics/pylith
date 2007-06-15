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
	tri3b.mesh \
	tri3c.mesh \
	tri3d.mesh \
	tri3_finalslip.spatialdb \
	tri3_sliptime.spatialdb \
	tri3_peakrate.spatialdb \
	tri3d_finalslip.spatialdb \
	tri3d_sliptime.spatialdb \
	tri3d_peakrate.spatialdb \
	quad4.mesh \
	quad4b.mesh \
	quad4c.mesh \
	quad4d.mesh \
	quad4e.mesh \
	quad4_finalslip.spatialdb \
	quad4_sliptime.spatialdb \
	quad4_peakrate.spatialdb \
	tet4.mesh \
	tet4b.mesh \
	tet4c.mesh \
	tet4d.mesh \
	tet4_finalslip.spatialdb \
	tet4_sliptime.spatialdb \
	tet4_peakrate.spatialdb \
	hex8.mesh \
	hex8b.mesh \
	hex8c.mesh \
	hex8d.mesh \
	hex8e.mesh \
	hex8f.mesh \
	hex8_finalslip.spatialdb \
	hex8_sliptime.spatialdb \
	hex8_peakrate.spatialdb

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
