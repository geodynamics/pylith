# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

TESTS = testcontrib.py
dist_check_SCRIPTS = testcontrib.py

noinst_PYTHON = \
	TestViscousFriction.py 


check-local: check-TESTS
	for f in $(TEST_LOGS); do printf "\n=== $$f ===\n\n"; cat $$f; done

# End of file 
