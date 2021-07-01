# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

dist_noinst_DATA = \
	slipfn.timedb \
	tri_a.mesh \
	tri_b.mesh \
	tri_c.mesh \
	tri_d.mesh \
	tri_e.mesh \
	tri_f.mesh \
	tri_g.mesh \
	tri_h.mesh \
	tri_i.mesh \
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
	quad_a.mesh \
	quad_b.mesh \
	quad_c.mesh \
	quad_d.mesh \
	quad_e.mesh \
	quad_f.mesh \
	quad_g.mesh \
	quad_h.mesh \
	quad_i.mesh \
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
	tet_a.mesh \
	tet_b.mesh \
	tet_c.mesh \
	tet_d.mesh \
	tet_e.mesh \
	tet_f.mesh \
	tet_g.mesh \
	tet_h.mesh \
	tet_i.mesh \
	tet_j.mesh \
	tet_k.jou \
	tet_k.mesh \
	tet4_finalslip.spatialdb \
	tet4_sliptime.spatialdb \
	tet4_risetime.spatialdb \
	tet4_sliprate.spatialdb \
	tet4e_finalslip.spatialdb \
	tet4e_sliptime.spatialdb \
	tet4e_risetime.spatialdb \
	tet4_initialtract.spatialdb \
	tet4_impulses.spatialdb \
	hex_a.mesh \
	hex_b.mesh \
	hex_c.mesh \
	hex_d.mesh \
	hex_e.mesh \
	hex_f.mesh \
	hex_g.mesh \
	hex_h.mesh \
	hex_i.mesh \
	hex_j.mesh \
	hex8_finalslip.spatialdb \
	hex8_sliptime.spatialdb \
	hex8_risetime.spatialdb \
	hex8_initialtract.spatialdb \
	hex8_impulses.spatialdb \
	exotoascii.py

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/libtests/faults/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data



# End of file 
