# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import make_suite
from pylith.topology.Field import Field
from pylith.topology.Mesh import Mesh


class TestField(unittest.TestCase):
    """Unit testing of Field object.
    """

    def test_constructor(self):
        mesh = Mesh()
        field = Field(mesh)
        self.assertTrue(not field is None)


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestField]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    unittest.main(verbosity=2)

    petsc.finalize()


# End of file
