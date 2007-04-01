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
##
## @brief Python PyLith abstract base class for solver.
##
## Factory: solver

from pyre.components.Component import Component

# Solver class
class Solver(Component):
  """
  Python abstract base class for solver.

  Factory: solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Solver facilities and properties.
    """

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
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solver")
    self.cppHandle = None
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return


# End of file 
