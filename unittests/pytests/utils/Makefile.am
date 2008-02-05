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

subpackage = utils
include $(top_srcdir)/subpackage.am

TESTS = testmanager.py testutils.py
check_SCRIPTS = testmanager.py testutils.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestEventLogger.py \
	TestPetscManager.py


# End of file 
