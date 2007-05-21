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

  Jacobian: A(t)
  solution: u(t+dt)
  residual: b(t) - A(t) \hat u(t+dt)
  constant: b(t)

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


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ExplicitElasticity import ExplicitElasticity
    return ExplicitElasticity()


  def initialize(self, mesh, materials, boundaryConditions, dimension, dt):
    """
    Initialize problem for explicit time integration.
    """
    self.integrators = []
    Formulation.initialize(self, mesh, materials, boundaryConditions,
                           dimension, dt)

    self._info.log("Initializing integrators.")

    self._info.log("Creating fields and matrices.")
    self.dispT = mesh.createRealSection("dispT", dimension)
    self.dispTmdt = mesh.createRealSection("dispTmdt", dimension)
    self.dispTpdt = mesh.createRealSection("dispTpdt", dimension)
    self.constant = mesh.createRealSection("constant", dimension)

    # Setup constraints
    # STUFF GOES HERE

    mesh.allocateRealSection(self.dispT)
    mesh.allocateRealSection(self.dispTmdt)
    mesh.allocateRealSection(self.dispTpdt)
    mesh.allocateRealSection(self.constant)
    
    self.jacobian = mesh.createMatrix(self.constant)

    self._info.log("Integrating Jacobian of operator.")
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, self.dispT)
    import pylith.utils.petsc as petsc
    petsc.mat_assemble(self.jacobian)

    self.solver.initialize(mesh, self.dispTpdt)
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    self._info.log("WARNING: Explicit::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second
    return dt
  

  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    self._info.log("WARNING: Explicit::prestep() not implemented.")
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("Integrating constant term in operator.")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(self.constant)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateConstant(self.constant, self.dispT, self.dispTmdt)

    self._info.log("Solving equations.")
    self.solver.solve(self.dispTpdt, self.jacobian, self.constant)
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    tmp = self.dispTmdt
    self.dispTmdt = self.dispT
    self.dispT = self.dispTpdt
    self.dispTpdt = tmp
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Explicit.
  """
  return Explicit()


# End of file 
