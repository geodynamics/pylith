#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#
# @file tests/pytests/topology/TestFields.py
#
# @brief Unit testing of Python Fields object.

import unittest

from pylith.topology.Fields import Fields
from pylith.topology.Mesh import Mesh


class TestFields(unittest.TestCase):
    """Unit testing of Fields object.
    """

    def test_constructor(self):
        mesh = Mesh()
        field = Fields(mesh)
        self.assertTrue(not field is None)


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFields))

    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    success = unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()

    petsc.finalize()


# End of file
