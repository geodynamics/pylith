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

## @file pylith/solver/Solver.py
## @brief Python PyLith abstract base class for solver.

from pyre.components.Component import Component

# Solver class
class Solver(Component):
  """Python abstract base class for solver."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Solver facilities and properties."""

    ## @class Inventory
    ## Python object for managing Solver facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solver"):
    """Constructor."""
    Component.__init__(self, name, facility="solver")
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    return


# version
__id__ = "$Id$"

# End of file 
