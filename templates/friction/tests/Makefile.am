# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
#

TESTS = testcontrib.py
dist_check_SCRIPTS = testcontrib.py

noinst_PYTHON = \
	TestViscousFriction.py 


check-local: check-TESTS
	for f in $(TEST_LOGS); do printf "\n=== $$f ===\n\n"; cat $$f; done

# End of file 
