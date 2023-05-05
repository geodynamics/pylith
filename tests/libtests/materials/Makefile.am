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

SUBDIRS = data

TESTS = test_materials

check_PROGRAMS = $(TESTS)


test_materials_SOURCES = \
	TestAuxiliaryFactoryElasticity.cc \
	TestAuxiliaryFactoryElasticity_Cases.cc \
	TestAuxiliaryFactoryElastic.cc \
	TestAuxiliaryFactoryElastic_Cases.cc \
	$(top_srcdir)/tests/src/FieldTester.cc \
	$(top_srcdir)/tests/src/driver_cppunit.cc


dist_noinst_HEADERS = \
	TestAuxiliaryFactoryElasticity.hh \
	TestAuxiliaryFactoryElastic.hh


# End of file
