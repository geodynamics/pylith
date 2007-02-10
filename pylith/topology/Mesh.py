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

## @file pylith/topology/Mesh.py
## @brief Python Mesh for finite-element topology information.

from pyre.components.Component import Component

# Mesh class
class Mesh(Component):
  """
  Python Mesh for finite-element topology information.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Mesh facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Mesh facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b coordsys Coordinate system associated with mesh

    import pyre.inventory

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pyre.inventory.facility("coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def distribute(self):
    """
    Distribute mesh across processors.
    """
    self._info.log("WARNING: Mesh::distribute() not implemented.")
    return self


  def __init__(self, name="mesh"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh")
    import pylith.topology.topology as bindings
    self.cppHandle = bindings.Mesh()
    return


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    self.coordsys = self.inventory.coordsys
    return
  

# End of file 
