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
	tri_small.mesh \
	tri_fourcells.mesh \
	quad_small.mesh \
	elasticstrain1d.spatialdb \
	tri3.mesh \
	tri3_disp.spatialdb \
	tri3_vel.spatialdb \
	tri3_disp2.spatialdb \
	tri3_vel2.spatialdb \
	tri3_tractions.spatialdb \
	tri3_force.spatialdb \
	tri3_force_rate.spatialdb \
	tri3_force_change.spatialdb \
	tri3_force.timedb \
	elasticplanestrain.spatialdb \
	quad4.mesh \
	quad4_disp.spatialdb \
	quad4_tractions.spatialdb \
	quad4_traction_initial.spatialdb \
	quad4_traction_rate.spatialdb \
	quad4_traction_change.spatialdb \
	quad4_traction.timedb \
	quad4_force.spatialdb \
	tet4.mesh \
	tet4_disp.spatialdb \
	tet4_disp2.spatialdb \
	tet4_vel2.spatialdb \
	tet4_tractions.spatialdb \
	tet4_force.spatialdb \
	elasticisotropic3d.spatialdb \
	hex8.mesh \
	hex8_disp.spatialdb \
	hex8_force.spatialdb \
	hex8b.mesh \
	hex8b_traction.spatialdb


# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/libtests/bc/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file
