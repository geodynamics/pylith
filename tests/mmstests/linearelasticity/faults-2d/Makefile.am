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

include $(top_srcdir)/tests/check_cppunit.am

TESTS = test_faultkin

check_PROGRAMS = $(TESTS)


test_faultkin_SOURCES = \
	TestFaultKin.cc \
	TestFaultKin2D_ConstRateDynamic.cc \
	$(top_srcdir)/tests/src/driver_cppunit.cc

#	TestFaultKin2D_TwoBlocksStatic.cc \
#	TestFaultKin2D_ThreeBlocksStatic.cc \
#	TestFaultKin2D_OneFaultShearNoSlip.cc \
#	TestFaultKin2D_TwoFaultsShearNoSlip.cc 


dist_noinst_HEADERS = \
	TestFaultKin.hh


dist_noinst_DATA = \
	data/tri.mesh \
	data/quad.mesh \
	data/tet.mesh \
	data/hex.mesh

noinst_TMP =

export_datadir = $(abs_builddir)


# End of file
