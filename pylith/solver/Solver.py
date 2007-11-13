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
    ## @li \b initial_guess_zero Use zero for initial guess.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    guessZero = pyre.inventory.bool("initial_guess_zero", default=True)
    guessZero.meta['tip'] = "Use zero for initial guess."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solver"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solver")
    self.cppHandle = None
    return


  def initialize(self, mesh, field):
    """
    Initialize solver.
    """
    assert(None != self.cppHandle)
    self.cppHandle.setInitialGuessNonzero(not self.guessZero)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.guessZero = self.inventory.guessZero
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solver():
  """
  Factory associated with Solver.
  """
  return Solver()


# End of file 
