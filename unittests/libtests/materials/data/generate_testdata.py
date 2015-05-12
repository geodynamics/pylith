#!/usr/bin/env python
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

## @file unittests/libtests/materials/data/generate_testdata.py

## @brief Python application for generating C++ data files for testing
## C++ objects.


tests = "all"

if tests in ["all", "elasticity", "elasticity-planestrain"]:
    import IsotropicLinearElasticityPlaneStrain; IsotropicLinearElasticityPlaneStrain.generate()

# End of file 
