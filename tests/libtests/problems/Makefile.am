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

TESTS = test_problems

check_PROGRAMS = test_problems

# Primary source files
test_problems_SOURCES = \
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
	$(top_srcdir)/tests/src/driver_cppunit.cc



dist_noinst_HEADERS = \
	TestPhysics.hh \
	TestObserversSoln.hh \
	TestObserversPhysics.hh \
	TestSolutionFactory.hh

noinst_TMP = \
	progress.txt \
	progress_time.txt \
	progress_step.txt

export_datadir = $(abs_builddir)


# End of file
