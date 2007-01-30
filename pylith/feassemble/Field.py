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

## @file pylith/feassemble/Field.py

## @brief Python PyLith field.

from pyre.components.Component import Component

# Field class
class Field(Component):
  """
  Python PyLith field.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Field facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Field facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """
    """
    return


  def __init__(self, name="field"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="field")
    self.sieveField = None
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////


  def _configure(self):
    """
    Set members based using inventory.
    """
    return
  

# End of file 
