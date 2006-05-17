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

## @file pylith/solver/SolverTSE.py
## @brief Python PyLith explicit time stepping solver.

from Solver import Solver

# SolverTSE class
class SolverTSE(Solver):
  """Python PyLith explicit time stepping solver."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Solver.Inventory):
    """Python object for managing SolverTSE facilities and properties."""

    ## @class Inventory
    ## Python object for managing SolverTSE facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solvertse"):
    """Constructor."""
    Solver.__init__(self, name)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    Solver._configure(self)
    return


# version
__id__ = "$Id$"

# End of file 
