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

SUBDIRS = data

TESTS = testbc.py

check_SCRIPTS = testbc.py

noinst_PYTHON = \
	TestDirichlet.py \
	TestBCSingle.py \
	TestBCTwoSides.py \
	TestBCFourSides.py \
	TestBCSixSides.py

TESTS_ENVIRONMENT = $(PYTHON)


# End of file 
