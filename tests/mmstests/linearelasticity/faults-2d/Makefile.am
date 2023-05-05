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
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/tests/check_catch2.am

TESTS = \
	run_twoblocksstatic2d.sh \
	run_threeblocksstatic2d.sh \
	run_onefaultshearnoslip2d.sh \
	run_twofaultsshearnoslip2d.sh

check_SCRIPTS = $(TESTS)

check_PROGRAMS = mmstest_linearelasticity_faults2d 


mmstest_linearelasticity_faults2d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestFaultKin.cc \
	TestCases.cc \
	TwoBlocksStatic.cc \
	ThreeBlocksStatic.cc \
	OneFaultShearNoSlip.cc \
	TwoFaultsShearNoSlip.cc


dist_noinst_HEADERS = \
	TestFaultKin.hh \
	TwoBlocksStatic.hh \
	ThreeBlocksStatic.hh \
	OneFaultShearNoSlip.hh \
	TwoFaultsShearNoSlip.hh


dist_noinst_DATA = \
	data/tri.mesh \
	data/quad.mesh \
	data/tet.mesh \
	data/hex.mesh


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
