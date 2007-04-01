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
##
## @brief Python abstract base class for formulations of solving equations.
##
## Factory: pde_formulation

from pyre.components.Component import Component

# Formulation class
class Formulation(Component):
  """
  Python abstract base class for formulations of solving equations.

  In general, we use some explicit or implicit formulation of the PDEs
  to create a linear form, [A]{u}={b} that we can solve.

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Formulation facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Formulation facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.

    import pyre.inventory

    from pylith.solver.SolverLinear import SolverLinear
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverLinear)
    solver.meta['tip'] = "Algebraic solver."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="formulation"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="pde_formulation")
    self.integrators = None
    return


  def initialize(self, mesh, materials, dimension, dt):
    """
    Create integrators for each element family.
    """
    raise NotImplementedError, "initialize() not implemented."
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.solver = self.inventory.solver
    return


  def _calcResidual(self):
    """
    Compute residual, {b(t)}.
    """
    raise NotImplementedError("Please implement calcResidual().")
    return


  def _calcJacobian(self):
    """
    Compute Jacobian, [A(t)].
    """
    raise NotImplementedError("Please implement calcJacobian().")
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Formulation.
  """
  return Formulation()


# End of file 
