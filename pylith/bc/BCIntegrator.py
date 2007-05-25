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

## @file pylith/bc/BCIntegrator.py
##
## @brief Python abstract base class for managing a boundary condition
## that requires integration.
##
## This implementation of a boundary condition applies to a single
## face of a domain and associates both a quadrature scheme with a
## physical boundary condition. Thus, applying different quadrature
## schemes along a face with the same physical boundary condition
## requires two boundary condition integrators, which can use the same
## database.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from pylith.feassemble.Integrator import Integrator

# BCIntegrator class
class BCIntegrator(BoundaryCondition, Integrator):
  """
  Python abstract base class for managing a boundary condition that
  requires integration.

  This implementation of a boundary condition applies to a single face
  of a domain and associates both a quadrature scheme with a physical
  boundary condition. Thus, applying different quadrature schemes
  along a face with the same physical boundary condition requires two
  boundary condition integrators, which can use the same database.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(BoundaryCondition.Inventory):
    """
    Python object for managing BCIntegrator facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BCIntegrator facilities and properties.
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

  def __init__(self, name="bcintegrator"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
    return


  def initialize(self, mesh):
    """
    Initialize boundary condition.
    """
    assert(None != self.cppHandle)
    Integrator.initQuadrature(self, self.quadrature)
    BoundaryCondition.initialize(self, mesh)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.quadrature = self.inventory.quadrature
    return

  
# End of file 
