# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

dist_noinst_DATA = \
	slipfn.timedb \
	static_friction_2d.spatialdb \
	static_friction_3d.spatialdb \
	tri3.mesh \
	tri3b.mesh \
	tri3c.mesh \
	tri3d.mesh \
	tri3e.mesh \
	tri3f.mesh \
	tri3g.mesh \
	tri3h.mesh \
	tri3_finalslip.spatialdb \
	tri3_finalslipB.spatialdb \
	tri3_sliptime.spatialdb \
	tri3_sliprate.spatialdb \
	tri3_risetime.spatialdb \
	tri3d_finalslip.spatialdb \
	tri3d_sliptime.spatialdb \
	tri3d_risetime.spatialdb \
	tri3g_finalslip.spatialdb \
	tri3g_sliptime.spatialdb \
	tri3g_risetime.spatialdb \
	tri3_initialtract.spatialdb \
	tri3_changetract.spatialdb \
	tri3d_initialtract.spatialdb \
	tri3_impulses.spatialdb \
	quad4.mesh \
	quad4b.mesh \
	quad4c.mesh \
	quad4d.mesh \
	quad4e.mesh \
	quad4f.mesh \
	quad4g.mesh \
	quad4h.mesh \
	quad4i.mesh \
	quad4_finalslip.spatialdb \
	quad4_sliptime.spatialdb \
	quad4_risetime.spatialdb \
	quad4e_finalslip.spatialdb \
	quad4e_sliptime.spatialdb \
	quad4e_risetime.spatialdb \
	quad4i_finalslip.spatialdb \
	quad4i_sliptime.spatialdb \
	quad4i_risetime.spatialdb \
	quad4_initialtract.spatialdb \
	quad4_impulses.spatialdb \
	tet4.mesh \
	tet4b.mesh \
	tet4c.mesh \
	tet4d.mesh \
	tet4e.mesh \
	tet4f.mesh \
	tet4g.mesh \
	tet4h.mesh \
	tet4i.mesh \
	tet4j.mesh \
	tet4k.jou \
	tet4k.mesh \
	tet4_finalslip.spatialdb \
	tet4_sliptime.spatialdb \
	tet4_risetime.spatialdb \
	tet4_sliprate.spatialdb \
	tet4e_finalslip.spatialdb \
	tet4e_sliptime.spatialdb \
	tet4e_risetime.spatialdb \
	tet4_initialtract.spatialdb \
	tet4_impulses.spatialdb \
	hex8.mesh \
	hex8b.mesh \
	hex8c.mesh \
	hex8d.mesh \
	hex8e.mesh \
	hex8f.mesh \
	hex8g.mesh \
	hex8h.mesh \
	hex8i.mesh \
	hex8j.mesh \
	hex8_finalslip.spatialdb \
	hex8_sliptime.spatialdb \
	hex8_risetime.spatialdb \
	hex8_initialtract.spatialdb \
	hex8_impulses.spatialdb \
	exotoascii.py

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/faults/data
export-data: $(dist_noinst_DATA)
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data



# End of file 
