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
# @file tests/pytests/mpi/TestReduce.py

# @brief Unit testing of MPI reduce functions.

import unittest

import pylith.mpi.mpi as mpi

# ----------------------------------------------------------------------
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


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestReduce))

    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    success = unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()

    petsc.finalize()


# End of file 
