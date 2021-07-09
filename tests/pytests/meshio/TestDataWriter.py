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
