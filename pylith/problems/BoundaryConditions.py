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

## @file pylith/problem/BoundaryConditions.py
##
## @brief Python container for boundary conditions.
##
## Factory: boundary_conditions

from pyre.components.Component import Component

# BoundaryConditions class
class BoundaryConditions(Component):
  """
  Python container for boundary conditions.

  Factory: boundary_conditions
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="boundaryconditions"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="boundary_conditions")
    self.bc = []
    return


# End of file 
