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

## @file pylith/problems/Explicit.py
##
## @brief Python Explicit object for solving equations using an
## explicit formulation.
##
## Factory: pde_formulation

from Formulation import Formulation

# Explicit class
class Explicit(Formulation):
  """
  Python Explicit object for solving equations using an explicit
  formulation.

  The formulation has the general form, [A(t)] {u(t+dt)} = {b(t)},
  where we want to solve for {u(t+dt)}, A(t) is usually constant
  (i.e., independent of time), and {b(t)} usually depends on {u(t)}
  and {u(t-dt)}.

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Formulation.Inventory):
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
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicit"):
    """
    Constructor.
    """
    Formulation.__init__(self, name)
    return


  def initialize(self, mesh, materialsBin):
    """
    Create explicit integrators for each element family.
    """
    self._info.log("Initializing integrators.")
    from pylith.feassemble.ExplicitElasticity import ExplicitElasticity
    
    self.integrators = []
    for material in materialsBin.materials:
      integrator = ExplicitElasticity()
      integrator.initQuadrature(material.quadrature)
      integrator.initMaterial(mesh, material)
      self.integrators.append(integrator)
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    self._info.log("WARNING: Explicit::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second
    return dt
  

  def prestep(self):
    """
    Hook for doing stuff before advancing time step.
    """
    self._info.log("WARNING: Explicit::prestep() not implemented.")
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("WARNING: Explicit::step() not implemented.")
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    self._info.log("WARNING: Explicit::poststep() not implemented.")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)
    return


  def _calcConstant(self):
    """
    Compute residual, {b(t)}.
    """
    self._info.log("WARNING: Explicit::calcConstant() not implemented.")
    return


  def _calcJacobian(self):
    """
    Compute Jacobian, [A(t)].
    """
    self._info.log("WARNING: Explicit::calcJacobian() not implemented.")
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Explicit.
  """
  return Explicit()


# End of file 
