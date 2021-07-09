#!/usr/bin/env python
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

## @file tests/pytests/mpi/TestReduce.py

## @brief Unit testing of MPI reduce functions.

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
