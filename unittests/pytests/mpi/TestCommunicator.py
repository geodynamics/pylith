#!/usr/bin/env python
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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/utils/TestCommunicator.py

## @brief Unit testing of Communicator object.

import unittest

import pylith.mpi.Communicator as mpicomm

# ----------------------------------------------------------------------
class TestCommunicator(unittest.TestCase):
  """
  Unit testing of Communicator object.
  """
  

  def test_petsc_comm_world(self):
    """
    Test petsc_comm_world().
    """
    comm = mpicomm.petsc_comm_world()
    return


  def test_petsc_comm_self(self):
    """
    Test petsc_comm_self().
    """
    comm = mpicomm.petsc_comm_self()
    return


  def test_mpi_comm_world(self):
    """
    Test mpi_comm_world().
    """
    comm = mpicomm.mpi_comm_world()
    return


  def test_mpi_comm_self(self):
    """
    Test mpi_comm_self().
    """
    comm = mpicomm.mpi_comm_self()
    return


  def test_rank(self):
    """
    Test Communicator.rank().
    """
    comm = mpicomm.petsc_comm_world()
    self.assertEqual(0, comm.rank)
    return


  def test_size(self):
    """
    Test Communicator.size().
    """
    comm = mpicomm.petsc_comm_world()
    self.assertEqual(1, comm.size)
    return


  def test_barrier(self):
    """
    Test Communicator.barrier().
    """
    comm = mpicomm.petsc_comm_world()
    comm.barrier()
    return


# End of file 
