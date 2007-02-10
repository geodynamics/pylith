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

## @file pylith/problems/Formulation.py

## @brief Python abstract base class for formulations of solving equations.

from pyre.components.Component import Component

# Formulation class
class Formulation(Component):
  """
  Python abstract base class for formulations of solving equations.

  In general, we use some explicit or implicit formulation of the PDEs
  to create a linear form, [A]{u}={b} that we can solve.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Explicit facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Explicit facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.

    import pyre.inventory

    #from pylith.solver.SolverTSE import SolverTSE
    #solver = pyre.inventory.facility("solver", factory=SolverTSE)
    #solver.meta['tip'] = "Algebraic solver."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """
    Initialize the formulation.
    """
    self._info.log("WARNING: Formulation::initialize not implemented.")
    return


  def calcResidual(self):
    """
    Compute residual, {b(t)}.
    """
    raise NotImplementedError("Please implement calcResidual().")
    return

  def calcJacobian(self):
    """
    Compute Jacobian, [A(t)].
    """
    raise NotImplementedError("Please implement calcJacobian().")
    return


  def __init__(self, name="formulation"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    #self.solver = self.inventory.solver
    return


# End of file 
