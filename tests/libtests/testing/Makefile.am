# -*- makefile -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


include $(top_srcdir)/tests/check_catch2.am

TESTS = libtest_testing

check_PROGRAMS = libtest_testing

# Primary source files
libtest_testing_SOURCES = \
	TestStubMethodTracker.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc

dist_noinst_HEADERS = 


# End of file
