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
    from pylith.feassemble.ElasticityExplicit import ElasticityExplicit
    return ElasticityExplicit()


  def initialize(self, mesh, materials, boundaryConditions,
                 interfaceConditions, dimension, dt):
    """
    Initialize problem for explicit time integration.
    """
    Formulation.initialize(self, mesh, materials, boundaryConditions,
                           interfaceConditions, dimension, dt)

    self._info.log("Creating fields and matrices.")
    self.fields.addReal("dispT")
    self.fields.addReal("dispTmdt")
    self.fields.addReal("dispTpdt")
    self.fields.addReal("residual")
    self.fields.createHistory(["dispTpdt", "dispT", "dispTmdt"])    
    self.fields.setFiberDimension("dispT", dimension)
    for constraint in self.constraints:
      constraint.setConstraintSizes(self.fields.getReal("dispT"), mesh)
    self.fields.allocate("dispT")
    for constraint in self.constraints:
      constraint.setConstraints(self.fields.getReal("dispT"), mesh)    
    self.fields.copyLayout("dispT")
    self.jacobian = mesh.createMatrix(self.fields.getReal("residual"))

    self._info.log("Forming Jacobian of operator.")
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(self.jacobian)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, self.fields)
    petsc.mat_assemble(self.jacobian)

    self.solver.initialize(mesh, self.fields.getReal("dispTpdt"))
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
    dispTpdt = self.fields.getReal("dispTpdt")
    for constraint in self.constraints:
      constraint.setField(t+dt, dispTpdt)

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
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("Integrating constant term in operator.")
    residual = self.fields.getReal("residual")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(residual)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateResidual(residual, self.fields)

    self._info.log("Solving equations.")
    self.solver.solve(self.fields.getReal("dispTpdt"), self.jacobian, residual)
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
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


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Explicit.
  """
  return Explicit()


# End of file 
