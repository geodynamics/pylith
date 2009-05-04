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

subpackage = bc
include $(top_srcdir)/subpackage.am

TESTS = testmpi.py
dist_check_SCRIPTS = testmpi.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestCommunicator.py


# End of file 
