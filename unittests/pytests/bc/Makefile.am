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
dist_check_SCRIPTS = testbc.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestAbsorbingDampers.py \
	TestDirichletBC.py \
	TestDirichletBoundary.py \
	TestNeumann.py


# End of file 
