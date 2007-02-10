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
    ## @li \b formulation Formulation for solving PDE

    import pyre.inventory

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facility("materials", factory=Homogeneous)
    materials.meta['tip'] = "Materials in problem."

    #from BoundaryConditions import BoundaryConditions
    #bc = pyre.inventory.facility("bc", factory=BoundaryConditions)
    #bc.meta['tip'] = "Boundary conditions."
  
    from Explicit import Explicit
    formulation = pyre.inventory.facility("formulation", factory=Explicit)
    formulation.meta['tip'] = "Formulation for solving PDE."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """
    Initialize problem by getting mesh, setting up boundary conditions, etc.
    """
    return


  def prestep(self):
    """
    User hook for doing stuff before advancing time step.
    """
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    return


  def poststep(self):
    """
    Update time and storage.
    """
    return


  def stableTimestep(self):
    """
    Determine stable time step for problem.
    """
    raise NotImplementedError, "Problem::stableTimestep() not implemented."
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    raise NotImplementedError, "Problem::checkpoint() not implemented."
    return
  

  def __init__(self, name="problem"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="problem")
    mesh = None
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
