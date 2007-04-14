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
	TestFault.py \
	TestFaultCohesive.py \
	TestFaultCohesiveKin.py


# End of file 
