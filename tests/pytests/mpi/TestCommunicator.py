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
import pylith.mpi.Communicator as mpicomm

class TestCommunicator(unittest.TestCase):
  """Unit testing of Communicator object.
  """
  

  def test_petsc_comm_world(self):
    comm = mpicomm.petsc_comm_world()
    return


  def test_petsc_comm_self(self):
    comm = mpicomm.petsc_comm_self()
    return


  def test_mpi_comm_world(self):
    comm = mpicomm.mpi_comm_world()
    return


  def test_mpi_comm_self(self):
    comm = mpicomm.mpi_comm_self()
    return


  def test_rank(self):
    comm = mpicomm.petsc_comm_world()
    self.assertEqual(0, comm.rank)
    return


  def test_size(self):
    comm = mpicomm.petsc_comm_world()
    self.assertEqual(1, comm.size)
    return


  def test_barrier(self):
    comm = mpicomm.petsc_comm_world()
    comm.barrier()
    return



def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestCommunicator]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    unittest.main(verbosity=2)

    petsc.finalize()


# End of file 
