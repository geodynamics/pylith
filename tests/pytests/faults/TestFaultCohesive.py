#!/usr/bin/env nemesis
# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


import unittest

from pylith.testing.TestCases import TestAbstractComponent, make_suite
from pylith.faults.FaultCohesive import FaultCohesive


class TestFaultCohesive(TestAbstractComponent):
    """Unit testing of FaultCohesive object."""

    _class = FaultCohesive


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestFaultCohesive]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
