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

SUBDIRS = data

TESTS = libtest_sources

check_PROGRAMS = $(TESTS)


libtest_sources_SOURCES = \
	TestAuxiliaryFactoryMomentTensorForce.cc \
	TestAuxiliaryFactoryMomentTensorForce_Cases.cc \
	TestAuxiliaryFactorySourceTime.cc \
	TestAuxiliaryFactorySourceTime_Cases.cc \
	TestAuxiliaryFactoryPointForce.cc \
	TestAuxiliaryFactoryPointForce_Cases.cc \
	$(top_srcdir)/tests/src/FieldTester.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


dist_noinst_HEADERS = \
	TestAuxiliaryFactoryMomentTensorForce.hh \
	TestAuxiliaryFactorySourceTime.hh \
	TestAuxiliaryFactoryPointForce.hh


# End of file
