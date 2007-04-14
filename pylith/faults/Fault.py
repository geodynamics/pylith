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

## @file pylith/faults/Fault.py
##

## @brief Python abstract base class for a fault surface.
##
## This implementation of a fault associates both physical
## properties and a quadrature scheme with the fault.
##
## Factory: fault

from pyre.components.Component import Component

# Fault class
class Fault(Component):
  """
  Python abstract base class for a fault surface.

  This implementation of a fault associates both physical
  properties and a quadrature scheme with the fault.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Fault facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Fault facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Fault identifier
    ## @li \b name Name of fault
    ##
    ## \b Facilities
    ## @li \b quadrature Quadrature object for numerical integration

    import pyre.inventory

    id = pyre.inventory.int("id", default=100)
    id.meta['tip'] = "Fault identifier (must be unique across all faults " \
                     "and materials)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name of material."

    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fault"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="fault")
    self.cppHandle = None
    return


  def adjustTopology(self, mesh):
    """
    Adjust mesh topology for fault implementation.
    """
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
    self.cppHandle.adjustTopology(mesh.cppHandle)
    return


  def initialize(self, mesh):
    """
    Initialize material property manager.
    """
    self._info.log("Initializing fault '%s'." % self.label)

    if 2 != self.quadrature.cell.cellDim:
      raise ValueError, \
            "Quadrature is incompatible with fault surface.\n" \
            "Dimensions for quadrature: %d, dimensions for surface: 2" % \
            self.quadrature.cell.cellDim

    #self.cppHandle.id = self.id
    #self.cppHandle.label = self.label
    #self.cppHandle.initialize(mesh.cppHandle, mesh.coordsys.cppHandle,
    #                          self.quadrature.cppHandle)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.id = self.inventory.id
    self.label = self.inventory.label
    self.quadrature = self.inventory.quadrature
    return

  
# End of file 
