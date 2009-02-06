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

## @file pylith/topology/Field.py
##
## @brief Python object for managing a vector field over vertices or
## cells of a finite-element mesh.

from topology import MeshField as ModuleMeshField
from topology import SubMeshField as ModuleSubMeshField

# ----------------------------------------------------------------------
# MeshField class
class MeshField(ModuleMeshField):
  """
  Python object for managing a vector field over vertices or cells of
  a finite-element mesh.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, mesh):
    """
    Constructor.
    """
    ModuleMeshField.__init__(self, mesh)
    return
    

# ----------------------------------------------------------------------
# SubMeshField class
class SubMeshField(ModuleSubMeshField):
  """
  Python object for managing a vector field over vertices or cells of
  a lower-dimension finite-element mesh.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, mesh):
    """
    Constructor.
    """
    ModuleSubMeshField.__init__(self, mesh)
    return
    

# End of file
