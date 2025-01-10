# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


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
