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

## @file pyre/problems/QuasiStatic.py
## @brief Python QuasiStatic application driver

from Problem import Problem

# QuasiStatic class
class QuasiStatic(Problem):
  """Python QuasiStatic for quasi-static crustal dynamics simulations."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Problem.Inventory):
    """Python object for managing QuasiStatic facilities and properties."""

    ## @class Inventory
    ## Python object for managing QuasiStatic facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.

    import pyre.inventory

    from pylith.solver.SolverTSI import SolverTSI
    solver = pyre.inventory.facility("solver", factory=SolverTSI)
    solver.meta['tip'] = "Algebraic solver."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    raise NotImplementedError, "QuasiStatic::initialize() not implemented."
    return


  def step(self):
    raise NotImplementedError, "QuasiStatic::step() not implemented."
    return


  def poststep(self):
    raise NotImplementedError, "QuasiStatic::poststep() not implemented."
    return


  def stableTimestep(self):
    raise NotImplementedError, "QuasiStatic::stableTimestep() not implemented."
    return


  def checkpoint(self):
    raise NotImplementedError, "QuasiStatic::checkpoint() not implemented."
    return
  

  def __init__(self, name="quasistatic"):
    """Constructor."""
    Problem.__init__(self, name)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

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
