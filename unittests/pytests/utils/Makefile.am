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

noinst_PYTHON = \
	TestEventLogger.py \
	TestPetscManager.py

TESTS_ENVIRONMENT = $(PYTHON)


# End of file 
