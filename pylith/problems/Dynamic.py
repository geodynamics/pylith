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

## @file pylith/problems/Dynamic.py
## @brief Python Dynamic for dynamic crustal dynamics simulations.

from Problem import Problem

# Dynamic class
class Dynamic(Problem):
  """Python Dynamic for dynamic crustal dynamics simulations."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Problem.Inventory):
    """Python object for managing Dynamic facilities and properties."""

    ## @class Inventory
    ## Python object for managing Dynamic facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.

    import pyre.inventory

    from pylith.solver.SolverTSE import SolverTSE
    solver = pyre.inventory.facility("solver", factory=SolverTSE)
    solver.meta['tip'] = "Algebraic solver."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    raise NotImplementedError, "Dynamic::initialize() not implemented."
    return


  def step(self):
    raise NotImplementedError, "Dynamic::step() not implemented."
    return


  def poststep(self):
    raise NotImplementedError, "Dynamic::poststep() not implemented."
    return


  def stableTimestep(self):
    raise NotImplementedError, "Dynamic::stableTimestep() not implemented."
    return


  def checkpoint(self):
    raise NotImplementedError, "Dynamic::checkpoint() not implemented."
    return
  

  def __init__(self, name="dynamic"):
    """Constructor."""
    Problem.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    self.solver = self.inventory.solver
    return


  def _calcResidual(self):
    """Compute solution residual."""
    return


  def _calcJacobian(self):
    """Calculation Jacobian."""
    return
  

# version
__id__ = "$Id$"

# End of file 
