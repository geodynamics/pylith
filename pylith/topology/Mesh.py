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
##
## Factory: finite_element_mesh

from topology import Mesh as ModuleMesh

# Mesh class
class Mesh(ModuleMesh):
  """
  Python Mesh for finite-element topology information.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, comm=None, dim=None):
    """
    Constructor.
    """
    if comm is None and dim is None:
      ModuleMesh.__init__(self)
    elif dim is None:
      ModuleMesh.__init__(self, comm.handle)
    else:
      ModuleMesh.__init__(self, comm.handle, dim)
    return


  def comm(self):
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
