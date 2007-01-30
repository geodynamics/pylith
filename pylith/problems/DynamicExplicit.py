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

## @file pylith/problems/DynamicExplicit.py

## @brief Python DynamicExplicit for dynamic simulations with explicit solver.

from Problem import Problem

# Dynamic class
class DynamicExplicit(Problem):
  """
  Python DynamicExplicit for dynamic simulations with explicit solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Problem.Inventory):
    """
    Python object for managing Dynamic facilities and properties.
    """

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
    """
    Create domain, bounday conditions, fields, and setup time loop.
    """
    Problem.initialize(self)
    
    #self.disptpdt.initialize()
    #self.dispt.initialize()
    #self.disptmdt.initialize()
    return


  def prestep(self):
    """
    User hook for doing stuff before advancing time step.
    """
    self._info.log("WARNING: DynamicExplicit::prestep not implemented.")
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("WARNING: DynamicExplicit::step not implemented.")
    return


  def poststep(self):
    """
    Update time and storage.
    """
    self._info.log("WARNING: DynamicExplicit::poststep not implemented.")
    return


  def stableTimestep(self):
    """
    Determine stable time step for problem.
    """
    self._info.log("WARNING: DynamicExplicit::stableTimestep not implemented.")
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    self._info.log("WARNING: DynamicExplicit::checkpoint not implemented.")
    return
  

  def __init__(self, name="dynamic"):
    """
    Constructor.
    """
    Problem.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Problem._configure(self)
    self.solver = self.inventory.solver
    self.disptpdt = self.inventory.disptpdt
    self.dispt = self.inventory.dispt
    self.disptmdt = self.inventory.disptmdt
    return


# End of file 
