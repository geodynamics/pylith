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
