# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = faults
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testfaults.py
dist_check_SCRIPTS = testfaults.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestBruneSlipFn.py \
	TestConstRateSlipFn.py \
	TestEqKinSrc.py \
	TestFaultCohesiveKin.py \
	TestFaultCohesiveDyn.py \
	TestTractPerturbation.py \
	TestFaultCohesiveImpulses.py \
	TestLiuCosSlipFn.py \
	TestSingleRupture.py \
	TestStepSlipFn.py \
	TestTimeHistorySlipFn.py


# End of file 
