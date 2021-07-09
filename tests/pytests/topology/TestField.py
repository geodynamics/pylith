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
# @file tests/pytests/topology/TestField.py
#
# @brief Unit testing of Python Field object.

import unittest

from pylith.topology.Field import Field
from pylith.topology.Mesh import Mesh


class TestField(unittest.TestCase):
    """Unit testing of Field object.
    """

    def test_constructor(self):
        mesh = Mesh()
        field = Field(mesh)
        self.assertTrue(not field is None)


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestField))

    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    success = unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()

    petsc.finalize()


# End of file
