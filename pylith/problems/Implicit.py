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
    from pylith.feassemble.ElasticityImplicit import ElasticityImplicit
    return ElasticityImplicit()


  def initialize(self, mesh, materials, boundaryConditions,
                 interfaceConditions, dimension, dt):
    """
    Initialize problem for implicit time integration.
    """
    Formulation.initialize(self, mesh, materials, boundaryConditions,
                           interfaceConditions, dimension, dt)

    self._info.log("Creating fields and matrices.")
    self.fields.addReal("dispTBctpdt")
    self.fields.addReal("dispIncr")
    self.fields.addReal("residual")
    self.fields.setFiberDimension("dispTBctpdt", dimension)
    for constraint in self.constraints:
      constraint.setConstraintSizes(self.fields.getReal("dispTBctpdt"))
    self.fields.allocate("dispTBctpdt")
    for constraint in self.constraints:
      constraint.setConstraints(self.fields.getReal("dispTBctpdt"))
    self.fields.copyLayout("dispTBctpdt")
    self.jacobian = mesh.createMatrix(self.fields.getReal("dispTBctpdt"))

    self.solver.initialize(mesh, self.fields.getReal("dispIncr"))

    from pyre.units.time import s
    self._solveElastic(mesh, materials, t=0.0*s, dt=dt)
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
    # Set dispTBctpdt to the BC at time t+dt. Unconstrained DOF are
    # unaffected and will be equal to their values at time t.
    dispTBctpdt = self.fields.getReal("dispTBctpdt")
    for constraint in self.constraints:
      constraint.setField(t+dt, dispTBctpdt)

    needNewJacobian = False
    for integrator in self.integrators:
      if integrator.needNewJacobian():
        needNewJacobian = True
    if needNewJacobian:
      self._info.log("Reforming Jacobian of operator.")
      import pylith.utils.petsc as petsc
      petsc.mat_setzero(self.jacobian)
      for integrator in self.integrators:
        integrator.timeStep(dt)
        integrator.integrateJacobian(self.jacobian, self.fields)
      petsc.mat_assemble(self.jacobian)
    # Put in loop over integrators to see if stiffness needs
    # reforming.  If so, then reform it.
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("Integrating residual term in operator.")
    residual = self.fields.getReal("residual")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(residual)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateResidual(residual, self.fields)

    self._info.log("Solving equations.")
    self.solver.solve(self.fields.getReal("dispIncr"), self.jacobian, residual)
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    # After solving, dispTBctpdt contains the displacements at time t
    # for unconstrained DOF and displacements at time t+dt at
    # constrained DOF. We add in the displacement increments (only
    # nonzero at unconstrained DOF) so that after poststrp(),
    # dispTBctpdt constains the solution at time t+dt.
    import pylith.topology.topology as bindings
    dispTBctpdt = self.fields.getReal("dispTBctpdt")
    bindings.addRealSections(dispTBctpdt, dispTBctpdt,
                             self.fields.getReal("dispIncr"))

    self._info.log("Updating integrators states.")
    for integrator in self.integrators:
      integrator.updateState(dispTBctpdt)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)
    return


  def _solveElastic(self, mesh, materials, t, dt):
    """
    Solve for elastic solution.
    """
    self._info.log("Computing elastic solution.")

    for material in materials.materials:
      material.useElasticBehavior(True)

    self._info.log("Setting constraints.")
    dispTBctpdt = self.fields.getReal("dispTBctpdt")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(dispTBctpdt)
    for constraint in self.constraints:
      constraint.setField(t, dispTBctpdt)

    self._info.log("Integrating Jacobian and residual of operator.")
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(self.jacobian)
    residual = self.fields.getReal("residual")
    dispIncr = self.fields.getReal("dispIncr")
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, self.fields)
      integrator.integrateResidual(residual, self.fields)
    import pylith.utils.petsc as petsc
    petsc.mat_assemble(self.jacobian)

    self._info.log("Solving equations.")
    self.solver.solve(dispIncr, self.jacobian, residual)

    import pylith.topology.topology as bindings
    dispTBctpdt = self.fields.getReal("dispTBctpdt")
    bindings.addRealSections(dispTBctpdt, dispTBctpdt, dispIncr)

    self._info.log("Updating integrators states.")
    for integrator in self.integrators:
      integrator.updateState(dispTBctpdt)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Implicit.
  """
  return Implicit()


# End of file 
