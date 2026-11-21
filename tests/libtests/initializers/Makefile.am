# -*- makefile -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


include $(top_srcdir)/tests/check_catch2.am

SUBDIRS = data

TESTS = libtest_initializers

check_PROGRAMS = libtest_initializers

libtest_initializers_SOURCES = \
	TestMeshReader.cc \
	TestMeshWriter.cc \
	TestMeshReordering.cc \
	TestMeshRefiner.cc \
	TestMeshDistributor.cc \
	TestMeshInsertInterfaces.cc \
	TestInitializer.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


clean-local:
	$(RM) $(RM_FLAGS) *.mesh


# End of file
