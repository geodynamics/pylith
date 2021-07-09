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
# @file tests/pytests/problems/TestSingleObserver.py
#
# @brief Unit testing of Python SingleObserver object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.problems.SingleObserver import (SingleSolnObserver, SinglePhysicsObserver)


class TestSingleSolnObserver(TestAbstractComponent):
    """Unit testing of SingleSolnObserver object.
    """
    _class = SingleSolnObserver


class TestSinglePhysicsObserver(TestAbstractComponent):
    """Unit testing of SinglePhysicsObserver object.
    """
    _class = SinglePhysicsObserver


if __name__ == "__main__":
    suite = unittest.TestSuite()
    classes = [
        TestSingleSolnObserver,
        TestSinglePhysicsObserver,
        ]
    for cls in classes:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
