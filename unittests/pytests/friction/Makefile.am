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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
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
	TestSlipWeakening.py \
	TestRateStateAgeing.py \
	TestTimeWeakening.py


# End of file 
