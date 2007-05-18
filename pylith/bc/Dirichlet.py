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

def validateDOF(value):
  """
  Validate list of fixed degrees of freedom.
  """
  try:
    size = len(value)
    num = map(float, value)
    for v in num:
      if v < 0:
        raise ValueError
  except:
    raise ValueError, \
          "'fixed_dof' must be a zero based list of indices of fixed " \
          "degrees of freedom."
  return value
  

# Dirichlet class
class Dirichlet(BoundaryCondition):
  """
  Python object for managing a Dirichlet (prescribed displacements)
  boundary condition.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(BoundaryCondition.Inventory):
    """
    Python object for managing BoundaryCondition facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BoundaryCondition facilities and properties.
    ##
    ## \b Properties
    ## @li \b fixed_dof Indices of fixed DOF (0=1st DOF, 1=2nd DOF, etc).
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    fixedDOF = pyre.inventory.list("fixed_dof", default=[],
                                   validator=validateDOF)
    fixedDOF.meta['tip'] = "Indices of fixed DOF (0=1st DOF, 1=2nd DOF, etc)."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dirichlet"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    self.fixedDOF = []
    import pylith.bc.bc as bindings
    self.cppHandle = bindings.Dirichlet()
    return


  def initialize(self, mesh):
    """
    Initialize Dirichlet boundary condition.
    """
    BoundaryCondition.initialize(self, mesh)
    self.cppHandle.fixedDOF = self.fixedDOF
    self.cppHandle.initialize(mesh.cppHandle, mesh.coordsys.cppHandle)
    return
  

  def setConstraintSizes(self, field, mesh):
    """
    Set number of constraints at points in field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.setConstraintSizes(field, mesh.cppHandle)
    return


  def setConstraints(self, field, mesh):
    """
    Set which degrees of freedom are constrained at points in field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.setConstraints(field, mesh.cppHandle)
    return


  def setField(self, t, field, mesh):
    """
    Set solution field at time t.
    """
    assert(None != self.cppHandle)
    self.cppHandle.setField(t, field, mesh.cppHandle)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.fixedDOF = self.inventory.fixedDOF
    return

  
# End of file 
