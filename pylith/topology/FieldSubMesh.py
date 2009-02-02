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

## @file pylith/topology/FieldSubMesh.py
##
## @brief Python object for managing a vector field over vertices or
## cells of a lower-dimension portion of a finite-element mesh.

from topology import FieldSubMesh as ModuleFieldSubMesh

# FieldSubMesh class
class Field(ModuleFieldSubMesh):
  """
  Python object for managing a vector field over vertices or cells of
  a lower-dimension portion of a finite-element mesh.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, submesh):
    """
    Constructor.
    """
    ModuleFieldSubMesh.__init__(self, submesh)
    return
    

# End of file
