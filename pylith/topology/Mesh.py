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

from topology import mpi_communicator

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
      print "A"
      ModuleMesh.__init__(self)
    elif dim is None:
      print "B"
      commInt = mpi_communicator(comm)
      print commInt
      ModuleMesh.__init__(self, comm)
    else:
      print "C"
      ModuleMesh.__init__(self, comm, dim)
    return


# End of file
