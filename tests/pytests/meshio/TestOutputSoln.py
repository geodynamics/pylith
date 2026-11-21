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
from pylith.meshio.OutputSoln import OutputSoln


class TestOutputSoln(TestAbstractComponent):
    """Unit testing of OutputSoln object."""

    _class = OutputSoln


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestOutputSoln]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
