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
    self.dispTBctpdt = mesh.createRealSection("dispTBctpdt", dimension)
    self.dispIncrement = mesh.createRealSection("dispIncrement", dimension)
    self.residual = mesh.createRealSection("residual", dimension)

    # Setup constraints
    # STUFF GOES HERE

    mesh.allocateRealSection(self.dispT)
    mesh.allocateRealSection(self.dispTBctpdt)
    mesh.allocateRealSection(self.dispIncrement)
    mesh.allocateRealSection(self.residual)

    self.jacobian = mesh.createMatrix(self.residual)

    self._info.log("Integrating Jacobian of operator.")
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, self.dispTBctpdt)
    import pylith.utils.petsc as petsc
    petsc.mat_assemble(self.jacobian)

    self.solver.initialize(mesh, self.dispIncrement)
    # self.mesh = mesh
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    self._info.log("WARNING: Implicit::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second
    return dt
  

  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    # This will need to set dispTBctpdt to the BC at time step t+dt.
    # Non-constrained DOF are unaffected and will be equal to their
    # values from time step t.
    # In this routine I also need to integrate the tractions for step
    # t+dt, but I don't think the function for this is available yet.
    # from pylith.bc.Dirichlet import Dirichlet

    # dispbc = Dirichlet

    # dispbc.setField(t+dt, self.dispTBctpdt, self.mesh)
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
      integrator.integrateResidual(self.residual, self.dispTBctpdt)

    self._info.log("Solving equations.")
    self.solver.solve(self.dispIncrement, self.jacobian, self.residual)
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    # This should give us the total displacements for time step t+dt, which
    # is renamed as time step t following the solve.
    # The vector dispTBctpdt contains the displacements from time step t
    # along with the displacement BC from time step t+dt.  The displacement
    # increments computed from the residual are then added to this to give us
    # the total displacement field at time t+dt.
    # Need a real way to do the operation below.
    # It is commented out for now, and should be replaced with a call to PETSc.
    # self.dispT = self.dispTBctpdt+self.dispIncrement
    self.dispTBctpdt = self.dispT
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
