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
    self.fields.addReal("dispT")
    self.fields.addReal("dispTBctpdt")
    self.fields.addReal("dispIncr")
    self.fields.addReal("residual")
    self.fields.createHistory(["dispTBctpdt", "dispT"])
    self.fields.setFiberDimension("dispT", dimension)
    for constraint in self.constraints:
      constraint.setConstraintSizes(self.fields.getReal("dispT"))
    self.fields.allocate("dispT")
    for constraint in self.constraints:
      constraint.setConstraints(self.fields.getReal("dispT"))
    self.fields.copyLayout("dispT")
    self.jacobian = mesh.createMatrix(self.fields.getReal("dispT"))

    from pyre.units.time import s
    self._solveElastic(mesh, t=0.0*s, dt=dt)

    self.solver.initialize(mesh, self.fields.getReal("dispIncr"))
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
        integrator.integrateJacobian(self.jacobian, self.fields.cppHandle)
      petsc.mat_assemble(self.jacobian)
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
      integrator.integrateResidual(residual, self.fields.cppHandle)

    self._info.log("Solving equations.")
    self.solver.solve(self.fields.getReal("dispIncr"), self.jacobian, residual)
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    # This should give us the total displacements after stepping
    # forward from t to t+dt, which is renamed as time step t for the
    # next solve. The field dispTBctpdt contains the displacements
    # from time step t along with the displacement BC from time step
    # t+dt.  The displacement increments computed from the residual
    # are then added to this to give us the total displacement field
    # at time t+dt.

    # Need a real way to do the operation below.
    # self.dispT = self.dispTBctpdt + self.dispIncr
    self.fields.shiftHistory()

    self._info.log("Updating integrators states.")
    for integrator in self.integrators:
      integrator.updateState(self.fields.getReal("dispT"))
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)
    return


  def _solveElastic(self, mesh, t, dt):
    """
    Solve for elastic solution.
    """
    self._info.log("Computing elastic solution.")

    self._info.log("Setting constraints.")
    dispT = self.fields.getReal("dispT")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(dispT)
    for constraint in self.constraints:
      constraint.setField(t, dispT)

    self._info.log("Integrating Jacobian and residual of operator.")
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(self.jacobian)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, self.fields.cppHandle)
      integrator.integrateResidual(self.fields.getReal("dispT"),
                                   self.fields.cppHandle)
    import pylith.utils.petsc as petsc
    petsc.mat_assemble(self.jacobian)

    self.solver.initialize(mesh, dispT)

    self._info.log("Solving equations.")
    self.solver.solve(dispT, self.jacobian, self.fields.getReal("residual"))

    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Implicit.
  """
  return Implicit()


# End of file 
