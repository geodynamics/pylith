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
from pylith.feassemble.Constraint import Constraint

def validateDOF(value):
  """
  Validate list of fixed degrees of freedom.
  """
  try:
    size = len(value)
    num = map(int, value)
    for v in num:
      if v < 0:
        raise ValueError
  except:
    raise ValueError, \
          "'fixed_dof' must be a zero based list of indices of fixed " \
          "degrees of freedom."
  return num
  

# Dirichlet class
class Dirichlet(BoundaryCondition, Constraint):
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
    Constraint.__init__(self)
    self._loggingPrefix = "DiBC "
    self.fixedDOF = []
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    Constraint.preinitialize(self, mesh)
    self.cppHandle.fixedDOF = self.fixedDOF    
    return


  def initialize(self):
    """
    Initialize Dirichlet boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    BoundaryCondition.initialize(self)

    self._logger.eventEnd(logEvent)    
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.fixedDOF = self.inventory.fixedDOF
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.bc.bc as bindings
      self.cppHandle = bindings.Dirichlet()    
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    Constraint._setupLogging(self)

    events = ["init"]
    for event in events:
      self._logger.registerEvent("%s%s" % (self._loggingPrefix, event))
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with Dirichlet.
  """
  return Dirichlet()

  
# End of file 
