# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


include $(top_srcdir)/tests/check_catch2.am

SUBDIRS = data

TESTS = libtest_problems

check_PROGRAMS = libtest_problems

# Primary source files
libtest_problems_SOURCES = \
	TestObserversSoln.cc \
	TestObserversPhysics.cc \
	TestSolutionFactory.cc \
	TestSolutionFactory_Cases.cc \
	TestProgressMonitor.cc \
	TestProgressMonitorTime.cc \
	TestProgressMonitorStep.cc \
	$(top_srcdir)/tests/src/ProgressMonitorStub.cc \
	$(top_srcdir)/tests/src/ObserverSolnStub.cc \
	$(top_srcdir)/tests/src/ObserverPhysicsStub.cc \
	$(top_srcdir)/tests/src/PhysicsImplementationStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/FieldTester.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc


dist_noinst_HEADERS = \
	TestPhysics.hh \
	TestSolutionFactory.hh

noinst_TMP = \
	progress.txt \
	progress_time.txt \
	progress_step.txt

export_datadir = $(abs_builddir)


# End of file
