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

subpackage = friction
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testfriction.py
dist_check_SCRIPTS = testfriction.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestFrictionModel.py \
	TestStaticFriction.py \
	TestSlipWeakening.py


# End of file 
