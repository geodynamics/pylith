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

TESTS = libtest_utils

check_PROGRAMS = libtest_utils

# Primary source files
libtest_utils_SOURCES = \
	TestEventLogger.cc \
	TestPyreComponent.cc \
	TestGenericComponent.cc \
	TestPylithVersion.cc \
	TestPetscVersion.cc \
	TestDependenciesVersion.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


dist_noinst_HEADERS = 


# End of file
