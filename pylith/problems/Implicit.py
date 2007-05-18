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

## @file pylith/problems/Implicit.py
##
## @brief Python Implicit object for solving equations using an
## implicit formulation.
##
## Factory: pde_formulation

from Formulation import Formulation

# Implicit class
class Implicit(Formulation):
  """
  Python Implicit object for solving equations using an implicit
  formulation.

  The formulation has the general form, [A(t+dt)] {u(t+dt)} = {b(t+dt)},
  where we want to solve for {u(t+dt)}. [A(t+dt)] generally
  depends on {u(t+dt)} as well as the current stresses and additional
  state variables.  For linear elastic or viscoelastic problems with
  constant time step size, A is a constant (after the elastic solution).
  {b(t+dt)} generally depends on the loads applied for time step t+dt
  (including the contributions to the internal force vector from
  displacement/velocity BC) as well as the internal force vector computed
  from the current stresses.

  Jacobian: A(t+dt)
  solution: u(t+dt)
  residual: bextern(t+dt) - bintern(t+dt)

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Formulation.Inventory):
    """
    Python object for managing Implicit facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Implicit facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="implicit"):
    """
    Constructor.
    """
    Formulation.__init__(self, name)
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ImplicitElasticity import ImplicitElasticity
    return ImplicitElasticity()


  def initialize(self, mesh, materials, boundaryConditions, dimension, dt):
    """
    Initialize problem for implicit time integration.
    """
    self.integrators = []
    Formulation.initialize(self, mesh, materials, boundaryConditions,
                           dimension, dt)

    self._info.log("Initializing integrators.")

    self._info.log("Creating fields and matrices.")
    self.dispT = mesh.createRealSection("dispT", dimension)
    self.dispTpdt = mesh.createRealSection("dispTpdt", dimension)
    self.residual = mesh.createRealSection("residual", dimension)

    # Setup constraints
    # STUFF GOES HERE

    mesh.allocateRealSection(self.dispT)
    mesh.allocateRealSection(self.dispTpdt)
    mesh.allocateRealSection(self.residual)

    self.jacobian = mesh.createMatrix(self.residual)

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
    self._info.log("WARNING: Implicit::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second
    return dt
  

  def prestep(self):
    """
    Hook for doing stuff before advancing time step.
    """
    self._info.log("WARNING: Implicit::prestep() not implemented.")
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("Integrating residual term in operator.")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(self.residual)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateResidual(self.residual, self.dispT)

    self._info.log("Solving equations.")
    self.solver.solve(self.dispTpdt, self.jacobian, self.residual)
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    self.dispT = self.dispTpdt
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
  Factory associated with Implicit.
  """
  return Implicit()


# End of file 
