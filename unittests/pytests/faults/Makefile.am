# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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
	TestLiuCosSlipFn.py \
	TestSingleRupture.py \
	TestStepSlipFn.py


# End of file 
