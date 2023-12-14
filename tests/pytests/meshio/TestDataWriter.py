#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/pytests/meshio/TestDataWriter.py
#
# @brief Unit testing of Python DataWriter object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.meshio.DataWriter import DataWriter


class TestDataWriter(TestAbstractComponent):
    """Unit testing of DataWriter object.
    """
    _class = DataWriter

    def test_mkfilename(self):
        writer = DataWriter()
        filename = writer.mkfilename(outputDir="abc", simName="defg", label="hijkl", suffix="hx3")
        self.assertEqual("abc/defg-hijkl.hx3", filename)


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDataWriter))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
