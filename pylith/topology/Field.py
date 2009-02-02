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

from topology import Field as ModuleField

# Field class
class Field(ModuleField):
  """
  Python object for managing a vector field over vertices or cells of
  a finite-element mesh.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, mesh):
    """
    Constructor.
    """
    ModuleField.__init__(self, mesh)
    return
    

# End of file
