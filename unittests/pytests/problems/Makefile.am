# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testproblems.py
dist_check_SCRIPTS = testproblems.py

noinst_PYTHON = \
	TestTimeStep.py \
	TestTimeStepAdapt.py \
	TestTimeStepUniform.py \
	TestTimeStepUser.py \
	TestProgressMonitor.py \
	TestProgressMonitorTime.py \
	TestProgressMonitorStep.py


# End of file 
