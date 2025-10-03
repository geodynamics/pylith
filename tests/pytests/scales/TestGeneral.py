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
from pylith.scales.General import General

from pythia.pyre.units.length import meter
from pythia.pyre.units.pressure import pascal
from pythia.pyre.units.time import second
from pythia.pyre.units.mass import kilogram
from pythia.pyre.units.temperature import kelvin


class TestGeneral(unittest.TestCase):

    def test_constructor(self):
        dim = General()
        dim._configure()

        self.assertEqual(1.0 * meter, dim.getDisplacementScale())
        self.assertEqual(1.0 * meter, dim.getLengthScale())
        self.assertEqual(1.0 * pascal, dim.getRigidityScale())
        self.assertEqual(1.0 * second, dim.getTimeScale())
        self.assertEqual(1.0 * kelvin, dim.getTemperatureScale())

    def test_displacementScale(self):
        dim = General()
        dim._configure()
        dim.setDisplacementScale(2.0 * meter)

        self.assertEqual(2.0 * meter, dim.getDisplacementScale())
        self.assertEqual(1.0 * pascal, dim.getLengthScale())
        self.assertEqual(1.0 * pascal, dim.getRigidityScale())
        self.assertEqual(1.0 * second, dim.getTimeScale())

    def test_lengthScale(self):
        dim = General()
        dim._configure()
        dim.setLengthScale(2.0 * meter)

        self.assertEqual(1.0 * meter, dim.getDisplacementScale())
        self.assertEqual(2.0 * meter, dim.getLengthScale())
        self.assertEqual(1.0 * pascal, dim.getRigidityScale())
        self.assertEqual(1.0 * second, dim.getTimeScale())

    def test_rigidityScale(self):
        dim = General()
        dim._configure()
        dim.setRigidityScale(2.0 * pascal)

        self.assertEqual(1.0 * meter, dim.getDisplacementScale())
        self.assertEqual(1.0 * meter, dim.getLengthScale())
        self.assertEqual(2.0 * pascal, dim.getRigidityScale())
        self.assertEqual(1.0 * second, dim.getTimeScale())

    def test_timeScale(self):
        dim = General()
        dim._configure()
        dim.setTimeScale(2.0 * second)

        self.assertEqual(1.0 * meter, dim.getDisplacementScale())
        self.assertEqual(1.0 * meter, dim.getLengthScale())
        self.assertEqual(1.0 * pascal, dim.getRigidityScale())
        self.assertEqual(2.0 * second, dim.getTimeScale())

    def test_temperatureScale(self):
        dim = General()
        dim._configure()
        dim.setTemperatureScale(2.0 * kelvin)

        self.assertEqual(1.0 * meter, dim.getDisplacementScale())
        self.assertEqual(1.0 * meter, dim.getLengthScale())
        self.assertEqual(1.0 * pascal, dim.getRigidityScale())
        self.assertEqual(1.0 * second, dim.getTimeScale())
        self.assertEqual(2.0 * kelvin, dim.getTemperatureScale())

    def test_nondimensionalize(self):
        dim = General()
        dim._configure()

        scale = 8.0 * meter
        value = 2.0 * meter
        valueE = 0.25

        self.assertEqual(valueE, dim.nondimensionalize(value, scale))

    def test_dimensionalize(self):
        dim = General()
        dim._configure()

        scale = 8.0 * meter
        value = 0.25
        valueE = 2.0 * meter

        self.assertEqual(valueE, dim.dimensionalize(value, scale))


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestGeneral]
    return make_suite(test_classes=TEST_CLASSES, loader=loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
