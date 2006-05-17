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
    ## @li \b disp_tpdt Displacement at time t+dt
    ## @li \b disp_t Displacement at time t
    ## @li \b disp_tmdt Displacement at time t-dt

    import pyre.inventory

    from pylith.solver.SolverTSE import SolverTSE
    solver = pyre.inventory.facility("solver", factory=SolverTSE)
    solver.meta['tip'] = "Algebraic solver."

    from pylith.feassemble.Field import Field
    disptpdt = pyre.inventory.facility("disp_tpdt", factory=Field,
                                       args=["disptpdt"])
    disptpdt.meta['tip'] = "Displacement at time t+dt."

    dispt = pyre.inventory.facility("disp_t", factory=Field,
                                       args=["dispt"])
    dispt.meta['tip'] = "Displacement at time t."

    disptmdt = pyre.inventory.facility("disp_tmdt", factory=Field,
                                       args=["disptmdt"])
    disptmdt.meta['tip'] = "Displacement at time t-dt."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    Problem.initialize(self)
    
    self.disptpdt.initialize()
    self.dispt.initialize()
    self.disptmdt.initialize()
    return


  def step(self, dt):
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
    Problem._configure(self)
    self.solver = self.inventory.solver
    self.disptpdt = self.inventory.disptpdt
    self.dispt = self.inventory.dispt
    self.disptmdt = self.inventory.disptmdt
    return


  def _calcResidual(self):
    """Compute solution residual."""
    return


  def _calcJacobian(self):
    """Calculation Jacobian."""
    return
  

  def _convergenceTest(self):
    """Test for convergence."""
    return

  
# version
__id__ = "$Id$"

# End of file 
