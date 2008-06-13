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
	TestFault.py \
	TestFaultCohesive.py \
	TestFaultCohesiveKin.py \
	TestEqKinSrc.py \
	TestStepSlipFn.py


# End of file 
