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

## @file pylith/bc/Neumann.py
##
## @brief Python object for managing traction boundary conditions.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from pylith.feassemble.Integrator import Integrator

# Neumann class
class Neumann(BoundaryCondition, Integrator):
  """
  Python object for managing traction boundary conditions.

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
    ## @li None
    ##
    ## \b Facilities
    ## @li \b quadrature Quadrature object for numerical integration

    import pyre.inventory

    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="neumann"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "NeBC "
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)
    self.quadrature.preinitialize()
    return


  def initialize(self):
    """
    Initialize Neumann boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix

    self._logger.eventBegin(logEvent)
    self.cppHandle.quadrature = self.quadrature.cppHandle
    BoundaryCondition.initialize(self)
    self._logger.eventEnd(logEvent)
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    BoundaryCondition.verifyConfiguration(self)
    if self.quadrature.cellDim != self.mesh.dimension()-1:
        raise ValueError, \
              "Quadrature scheme and mesh are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension of mesh boundary '%s': %d" % \
              (self.quadrature.cellDim,
               self.label, self.mesh.dimension()-1)    
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.quadrature = self.inventory.quadrature
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.bc.bc as bindings
      self.cppHandle = bindings.Neumann()    
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    Integrator._setupLogging(self)

    events = ["init"]
    for event in events:
      self._logger.registerEvent("%s%s" % (self._loggingPrefix, event))
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with Neumann.
  """
  return Neumann()

  
# End of file 
