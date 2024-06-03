# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import make_suite
import pylith.mpi.mpi as mpi

class TestReduce(unittest.TestCase):
  """Unit testing of MPI reduce functions.
  """
  
  def test_allreduce_scalar_double(self):
    value = 2.0
    result = mpi.allreduce_scalar_double(value, mpi.mpi_sum(), mpi.petsc_comm_world())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_double(value, mpi.mpi_min(), mpi.petsc_comm_self())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_double(value, mpi.mpi_max(), mpi.petsc_comm_world())
    self.assertEqual(value, result)
    return


  def test_allreduce_scalar_int(self):
    value = 3
    result = mpi.allreduce_scalar_int(value, mpi.mpi_sum(), mpi.petsc_comm_world())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_int(value, mpi.mpi_min(), mpi.petsc_comm_self())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_int(value, mpi.mpi_max(), mpi.petsc_comm_world())
    self.assertEqual(value, result)
    return


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestReduce]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    unittest.main(verbosity=2)

    petsc.finalize()


# End of file 
