#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pylith/mpi/Communicator.py
##
## @brief Python MPI communicator object.
##
## Provides SWIG friendly interface to MPI communicator object.

import pylith.mpi.mpi as mpimodule

# Communicator class
class Communicator(object):
  """
  Python MPI communicator object.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, cppComm):
    """
    Constructor.
    """
    self.handle = cppComm
    self.rank = mpimodule.rank(self.handle)
    self.size = mpimodule.size(self.handle)
    return


  def barrier(self):
    """
    MPI Barrier.
    """
    mpimodule.barrier(self.handle)
    return


# ----------------------------------------------------------------------
def petsc_comm_world():
  """
  Python wrapper around PETSC_COMM_WORLD.
  """
  global _petsc_world
  if _petsc_world is None:
      _petsc_world = Communicator(mpimodule.petsc_comm_world())
  return _petsc_world


# ----------------------------------------------------------------------
def petsc_comm_self():
  """
  Python wrapper around PETSC_COMM_SELF.
  """
  global _petsc_self
  if _petsc_self is None:
      _petsc_self = Communicator(mpimodule.petsc_comm_self())
  return _petsc_self


# ----------------------------------------------------------------------
def mpi_comm_world():
  """
  Python wrapper around MPI_COMM_WORLD.
  """
  global _mpi_world
  if _mpi_world is None:
      _mpi_world = Communicator(mpimodule.mpi_comm_world())
  return _mpi_world


# ----------------------------------------------------------------------
def mpi_comm_self():
  """
  Python wrapper around MPI_COMM_SELF.
  """
  global _mpi_self
  if _mpi_self is None:
      _mpi_self = Communicator(mpimodule.mpi_comm_self())
  return _mpi_self


# ----------------------------------------------------------------------
# Singletons
_petsc_world = None
_petsc_self = None
_mpi_world = None
_mpi_self = None


# End of file
