#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/topology/Mesh.py
##
## @brief Python Mesh for finite-element topology information.

from topology import Mesh as ModuleMesh

# Mesh class
class Mesh(ModuleMesh):
  """
  Python Mesh for finite-element topology information.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, dim=None, comm=None):
    """
    Constructor.
    """
    if comm is None and dim is None:
      ModuleMesh.__init__(self)
    elif comm is None:
      ModuleMesh.__init__(self, dim)
    else:
      ModuleMesh.__init__(self, dim, comm.handle)

    # Name of logging stage for mesh. We progress through various
    # stages as we read, distribute, and refine mesh.
    self.memLoggingStage = "Mesh"
    return


  def setComm(self, comm):
    """
    Set communicator.
    """
    ModuleMesh.comm(self, comm.handle)
    return


  def getComm(self):
    """
    Get communicator.
    """
    # Use Communicator object to wrap C++ MPI_Comm* returned by
    # module.
    from pylith.mpi.Communicator import Communicator
    return Communicator(ModuleMesh.comm(self))


  def checkMaterialIds(self, materialIds):
    """
    Check material ids for consistency with mesh.
    """
    from topology import MeshOps_checkMaterialIds
    MeshOps_checkMaterialIds(self, materialIds)
    return


  def dimension(self):
    """
    Return the topological mesh dimension
    """
    return ModuleMesh.dimension(self)


  def coneSize(self):
    """
    Return the representative cone size, or number of vertices on a cell
    """
    return ModuleMesh.coneSize(self)


  def numVertices(self):
    """
    Return the number of vertices.
    """
    return ModuleMesh.numVertices(self)


  def numCells(self):
    """
    Return the number of cells.
    """
    return ModuleMesh.numCells(self)

  def groupSizes(self):
    """
    Return the name and number of vertices for each group
    """
    groups = []
    names  = ModuleMesh.groups(self)
    for name in names:
      groups.append((name,ModuleMesh.groupSize(self, name)))
    return groups

# End of file
