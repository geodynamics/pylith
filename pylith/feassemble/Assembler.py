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

## @file pyre/feassemble/Assembler.py
## @brief Python finite-element assembler.

from pyre.components.Component import Component

# Assembler class
class Assembler(Component):
  """Python finite-element assembler."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Assembler facilities and properties."""

    ## @class Inventory
    ## Python object for managing Assembler facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="assembler"):
    """Constructor."""
    Component.__init__(self, name)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    return
  

# version
__id__ = "$Id$"

# End of file 
