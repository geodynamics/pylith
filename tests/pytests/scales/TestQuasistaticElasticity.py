#!/usr/bin/env nemesis
#
# =================================================================================================
# This code is part of SpatialData, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/spatialdata).
#
# Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

import unittest

from spatialdata.testing.TestCases import make_suite
from pylith.scales.QuasistaticElasticity import QuasistaticElasticity

from pythia.pyre.units.length import meter, kilometer
from pythia.pyre.units.pressure import pascal
from pythia.pyre.units.time import year


class TestElasticityScales(unittest.TestCase):

    def test_constructor(self):
        dim = QuasistaticElasticity()
        dim._configure()

        # Default values
        displacementScale = 1.0 * meter
        lengthScale = 100.0 * kilometer
        rigidityScale = 1.0e10 * pascal
        timeScale = 1.0e2 * year

        # Check defaults
        self.assertEqual(displacementScale, dim.getDisplacementScale())
        self.assertEqual(lengthScale, dim.getLengthScale())
        self.assertEqual(rigidityScale, dim.getRigidityScale())
        self.assertEqual(timeScale, dim.getTimeScale())


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestElasticityScales]
    return make_suite(test_classes=TEST_CLASSES, loader=loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
