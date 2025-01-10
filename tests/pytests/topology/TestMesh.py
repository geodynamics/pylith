# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import make_suite
from pylith.topology.Mesh import Mesh


class TestMesh(unittest.TestCase):
    """Unit testing of Mesh object.
    """

    def test_constructor(self):
        mesh = Mesh()
        self.assertTrue(not mesh is None)


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestMesh]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    unittest.main(verbosity=2)

    petsc.finalize()


# End of file
