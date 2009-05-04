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

from pylith.utils.PetscComponent import PetscComponent

# Solver class
class Solver(PetscComponent):
  """
  Python abstract base class for solver.

  Factory: solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
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

    guessZero = pyre.inventory.bool("initial_guess_zero", default=True)
    guessZero.meta['tip'] = "Use zero for initial guess."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solver"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="solver")
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solver():
  """
  Factory associated with Solver.
  """
  return Solver()


# End of file 
