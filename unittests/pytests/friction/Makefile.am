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

TESTS = testfriction.py
dist_check_SCRIPTS = testfriction.py

noinst_PYTHON = \
	TestFrictionModel.py \
	TestStaticFriction.py \
	TestSlipWeakening.py \
	TestSlipWeakeningTime.py \
	TestSlipWeakeningTimeStable.py \
	TestRateStateAgeing.py \
	TestTimeWeakening.py


# End of file 
