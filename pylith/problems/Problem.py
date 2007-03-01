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

## @file pylith/problems/Problem.py

## @brief Python abstract base class for crustal dynamics problems.

from pyre.components.Component import Component

# Problem class
class Problem(Component):
  """
  Python abstract base class for crustal dynamics problems.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Problem facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Problem facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b materials Materials in problem.
    ## @li \b bc Boundary conditions.

    import pyre.inventory

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facility("materials", factory=Homogeneous)
    materials.meta['tip'] = "Materials in problem."

    #from BoundaryConditions import BoundaryConditions
    #bc = pyre.inventory.facility("bc", factory=BoundaryConditions)
    #bc.meta['tip'] = "Boundary conditions."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="problem"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="problem")
    self.mesh = None
    return


  def initialize(self, mesh):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    raise NotImplementedError, "initialize() not implemented."
    return


  def run(self, app):
    """
    Solve the problem.
    """
    raise NotImplementedError, "run() not implemented."
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    raise NotImplementedError, "checkpoint() not implemented."
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    self.materials = self.inventory.materials
    #self.bc = self.inventory.bc
    self.formulation = self.inventory.formulation
    return


# End of file 
