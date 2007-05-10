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

## @file pylith/bc/Dirichlet.py
##
## @brief Python object for managing a Dirichlet (prescribed
## displacements) boundary condition.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition

# Dirichlet class
class Dirichlet(BoundaryCondition):
  """
  Python object for managing a Dirichlet (prescribed displacements)
  boundary condition.

  Factory: boundary_condition
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dirichlet"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    #import pylith.bc.bc as bindings
    #self.cppHandle = bindings.Dirichlet()
    return


  def initialize(self, mesh):
    """
    Initialize Dirichlet boundary condition.
    """
    return
  

  def setField(self, field, t):
    """
    Set solution field at time t.
    """
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    return

  
# End of file 
