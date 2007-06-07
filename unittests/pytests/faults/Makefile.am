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

check_SCRIPTS = testfaults.py

noinst_PYTHON = \
	TestBruneSlipFn.py \
	TestFault.py \
	TestFaultCohesive.py \
	TestFaultCohesiveKin.py

TESTS_ENVIRONMENT = $(PYTHON)


# End of file 
