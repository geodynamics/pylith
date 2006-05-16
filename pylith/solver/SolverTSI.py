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

## @file pylith/solver/SolverTSI.py
## @brief Python implicit time stepping solver.

from pyre.components.Component import Component

# SolverTSI class
class SolverTSI(Component):
  """Python implicit time stepping solver."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Problem.Inventory):
    """Python object for managing SolverTSI facilities and properties."""

    ## @class Inventory
    ## Python object for managing SolverTSI facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quasistatic"):
    """Constructor."""
    Problem.__init__(self, name)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    return

# version
__id__ = "$Id$"

# End of file 
