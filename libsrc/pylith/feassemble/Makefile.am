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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = feassemble
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	IntegratorPointwise.hh \
	ConstraintPointwise.hh \
	Observer.hh \
	ObservedComponent.hh \
	AuxiliaryFactory.hh \
	feassemblefwd.hh

noinst_HEADERS =

dist_noinst_DATA = \
	jacobian2d_lgdeform.wxm \
	jacobian2d_nonsymm_lgdeform.wxm \
	jacobian3d_lgdeform.wxm \
	jacobian3d_nonsymm_lgdeform.wxm \
	tri3_elasticity.wxm \
	tet4_elasticity.wxm


# End of file
