# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


import unittest

from pylith.testing.TestCases import TestComponent, make_suite
from pylith.materials.Elasticity import Elasticity, material


class TestElasticity(TestComponent):
    """Unit testing of Elasticity object."""

    _class = Elasticity
    _factory = material


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestElasticity]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
