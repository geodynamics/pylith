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

## @file pylith/bc/AbsorbingDampers.py
##
## @brief Python object for managing absorbing boundary condition
## using simple dashpots.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from pylith.feassemble.Integrator import Integrator

# AbsorbingDampers class
class AbsorbingDampers(BoundaryCondition, Integrator):
  """
  Python object for managing absorbing boundary condition using simple
  dashpots.

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

  def __init__(self, name="absorbingdampers"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    self.quadrature.preinitialize()
    return


  def initialize(self):
    """
    Initialize AbsorbingDampers boundary condition.
    """
    self.cppHandle.quadrature = self.quadrature.cppHandle
    BoundaryCondition.initialize(self)
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    BoundaryCondition.verifyConfiguration(self)
    if self.quadrature.spaceDim != self.mesh.dimension-1:
        raise ValueError, \
              "Quadrature scheme and mesh are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension of mesh boundary '%s': %d" % \
              (self.quadrature.spaceDim, self.label, self.mesh.dimension-1)    
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.bc.bc as bindings
      self.cppHandle = bindings.AbsorbingDampers()    
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with AbsorbingDampers.
  """
  return AbsorbingDampers()

  
# End of file 
