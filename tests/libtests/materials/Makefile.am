# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


include $(top_srcdir)/tests/check_catch2.am

SUBDIRS = data

TESTS = libtest_materials

check_PROGRAMS = $(TESTS)


libtest_materials_SOURCES = \
	TestAuxiliaryFactoryElasticity.cc \
	TestAuxiliaryFactoryElasticity_Cases.cc \
	TestAuxiliaryFactoryLinearElastic.cc \
	TestAuxiliaryFactoryLinearElastic_Cases.cc \
	$(top_srcdir)/tests/src/FieldTester.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


dist_noinst_HEADERS = \
	TestAuxiliaryFactoryElasticity.hh \
	TestAuxiliaryFactoryLinearElastic.hh


# End of file
