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
  

# End of file
