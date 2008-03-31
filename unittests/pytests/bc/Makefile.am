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

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestAbsorbingDampers.py \
	TestNeumann.py \
	TestDirichletPoints.py \
	TestDirichletBoundary.py


# End of file 
