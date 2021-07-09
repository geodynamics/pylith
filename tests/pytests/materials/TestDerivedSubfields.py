#!/usr/bin/env nemesis
#
# ======================================================================
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
# ======================================================================
#
# @file tests/pytests/materials/TestDerivedSubfields.py
#
# @brief Unit testing of Python derived subfields object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
import pylith.materials.DerivedSubfieldsElasticity


class TestDerivedSubfieldsElasticity(TestComponent):
    """Unit testing of DerivedSubfieldsElasticity object.
    """
    _class = pylith.materials.DerivedSubfieldsElasticity.DerivedSubfieldsElasticity
    _factory = pylith.materials.DerivedSubfieldsElasticity.derived_subfields


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestDerivedSubfieldsElasticity,
        ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
