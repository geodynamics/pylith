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

## @file unittests/pytests/utils/TestReduce.py

## @brief Unit testing of MPI reduce functions.

import unittest

import pylith.mpi.mpi as mpi

# ----------------------------------------------------------------------
class TestReduce(unittest.TestCase):
  """
  Unit testing of MPI reduce functions.
  """
  

  def test_allreduce_scalar_double(self):
    """
    Test allreduce_double().
    """
    value = 2.0
    result = mpi.allreduce_scalar_double(value, mpi.mpi_sum(), mpi.petsc_comm_world())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_double(value, mpi.mpi_min(), mpi.petsc_comm_self())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_double(value, mpi.mpi_max(), mpi.petsc_comm_world())
    self.assertEqual(value, result)
    return


  def test_allreduce_scalar_int(self):
    """
    Test allreduce_int().
    """
    value = 3
    result = mpi.allreduce_scalar_int(value, mpi.mpi_sum(), mpi.petsc_comm_world())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_int(value, mpi.mpi_min(), mpi.petsc_comm_self())
    self.assertEqual(value, result)

    result = mpi.allreduce_scalar_int(value, mpi.mpi_max(), mpi.petsc_comm_world())
    self.assertEqual(value, result)
    return


# End of file 
