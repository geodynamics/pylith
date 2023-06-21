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

TESTS = libtest_bc

check_PROGRAMS = libtest_bc

# Primary source files
libtest_bc_SOURCES = \
	Placeholder.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc

#	TestDirichletTimeDependent.cc \
#	TestDirichletTimeDependent_Cases.cc \
#	TestNeumannTimeDependent.cc \
#	TestNeumannTimeDependent_Cases.cc \
#	TestAbsorbingDampers.cc \
#	TestAbsorbingDampers_Cases.cc


dist_noinst_HEADERS = \
	TestDirichletTimeDependent.hh \
	TestNeumannTimeDependent.hh \
	TestAbsorbingDampers.hh


# End of file
