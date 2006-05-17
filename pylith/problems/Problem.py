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
  """Python abstract base class for crustal dynamics problems."""
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Problem facilities and properties."""

    ## @class Inventory
    ## Python object for managing Problem facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b mesh Finite-element topology.
    ## @li \b assembler Finite-element assembler.
    ## @li \b coordinates Field for coordinates of vertices associated
    ##        with domain
    ## @li \b equation_type Field defining different constitutive models
    ## @li \b materials Materials in problem.

    import pyre.inventory

    from pylith.topology.Mesh import Mesh
    mesh = pyre.inventory.facility("mesh", factory=Mesh)
    mesh.meta['tip'] = "Finite-element topology."

    from pylith.feassemble.Assembler import Assembler
    assembler = pyre.inventory.facility("assembler", factory=Assembler)
    assembler.meta['tip'] = "Finite-element assembler."

    from pylith.feassemble.Field import Field
    coordinates = pyre.inventory.facility("coordinates", factory=Field,
                                          args=["coordinates"])
    coordinates.meta['tip'] = "Field for coordinates of vertices associated" \
                              "with domain"

    eqntype = pyre.inventory.facility("equation_type", factory=Field,
                                      args=["eqntype"])
    eqntype.meta['tip'] = "Field defining different constitutive models."

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facility("materials", factory=Homogeneous)
    mesh.meta['tip'] = "Materials in problem."
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """Create domain, bounday conditions, fields, and setup time loop."""
    self.coordinates.initialize()
    self.eqntype.initialize()
    return


  def prestep(self):
    """User hook for doing stuff before advancing time step."""
    return


  def step(self, dt):
    """Advance to next time step."""
    return


  def poststep(self):
    """Update time and storage."""
    return


  def stableTimestep(self):
    """Determine stable time step for problem."""
    raise NotImplementedError, "Problem::stableTimestep() not implemented."
    return


  def checkpoint(self):
    """Save problem state for restart."""
    raise NotImplementedError, "Problem::checkpoint() not implemented."
    return
  

  def __init__(self, name="problem"):
    """Constructor."""
    Component.__init__(self, name, facility="problem")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    self.mesh = self.inventory.mesh
    self.assembler = self.inventory.assembler
    self.coordinates = self.inventory.coordinates
    self.eqntype = self.inventory.eqntype
    return

# version
__id__ = "$Id$"

# End of file 
