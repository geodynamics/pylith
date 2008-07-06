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

subpackage = problems
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testproblems.py
dist_check_SCRIPTS = testproblems.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestTimeStepUniform.py


# End of file 
