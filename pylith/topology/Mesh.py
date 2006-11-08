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
  """Python Mesh for finite-element topology information."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Mesh facilities and properties."""

    ## @class Inventory
    ## Python object for managing Mesh facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory


  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def distribute(self):
    #raise NotImplementedError, "Mesh::distribute() not implemented."
    return self


  def __init__(self, name="mesh"):
    """Constructor."""
    Component.__init__(self, name, facility="mesh")
    self.cppHandle = None
    return


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    return
  

# End of file 
