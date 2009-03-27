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
from pylith.utils.profiling import resourceUsageString

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
    self._loggingPrefix = "TSIm "
    self.solnField = {'name': "disp(t), bc(t+dt)",
                      'label': "displacement"}
    self._stepCount = None
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityImplicit import ElasticityImplicit
    return ElasticityImplicit()


  def initialize(self, dimension, normalizer):
    """
    Initialize problem for implicit time integration.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    Formulation.initialize(self, dimension, normalizer)

    self._info.log("Creating other fields.")
    self.fields.add("dispIncr(t)", "displacement increment")
    self.fields.add("residual", "residual")
    self.fields.copyLayout("disp(t), bc(t+dt)")
    self._debug.log(resourceUsageString())

    self._info.log("Creating Jacobian matrix.")
    from pylith.topology.Jacobian import Jacobian
    self.jacobian = Jacobian(self.fields)
    self.jacobian.zero() # TEMPORARY, to get correct memory usage
    self._debug.log(resourceUsageString())

    # Create Petsc vectors for fields involved in solve
    dispIncr = self.fields.get("dispIncr(t)")
    dispIncr.createVector()
    residual = self.fields.get("residual")
    residual.createVector()

    self._info.log("Initializing solver.")
    self.solver.initialize(self.fields, self.jacobian, self)
    self._debug.log(resourceUsageString())

    # Initial time step solves for total displacement field, not increment
    self._stepCount = 0
    for constraint in self.constraints:
      constraint.useSolnIncr(False)
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.useSolnIncr(False)

    self._logger.eventEnd(logEvent)
    return


  def getStartTime(self):
    """
    Get time at which time stepping should start.
    """
    dt = self.timeStep.currentStep()
    return -dt


  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    logEvent = "%sprestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    # Set dispTBctpdt to the BC t time t+dt. Unconstrained DOF are
    # unaffected and will be equal to their values at time t.
    self._info.log("Setting constraints.")
    dispTBctpdt = self.fields.get("disp(t), bc(t+dt)")
    for constraint in self.constraints:
      constraint.setField(t+dt, dispTBctpdt)

    # If finishing first time step, then switch from solving for total
    # displacements to solving for incremental displacements
    if 1 == self._stepCount:
      self._info.log("Switching from total field solution to incremental " \
                     "field solution.")
      for constraint in self.constraints:
        constraint.useSolnIncr(True)
      for integrator in self.integratorsMesh + self.integratorsSubMesh:
        integrator.useSolnIncr(True)
      self._reformJacobian(t, dt)

    ### NONLINEAR: Might want to move logic into IntegrateJacobian() and set a flag instead
    needNewJacobian = False
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.timeStep(dt)
      if integrator.needNewJacobian():
        needNewJacobian = True
    if needNewJacobian:
      self._reformJacobian(t, dt)

    self._logger.eventEnd(logEvent)
    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    logEvent = "%sstep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    dispIncr = self.fields.get("dispIncr(t)")
    dispIncr.zero()

    ### NONLINEAR: This moves under SNES control as IntegrateResidual()
    ### NONLINEAR: Also move updateState() from Integrator.poststep() to this function
    self._reformResidual(t+dt, dt)

    self._info.log("Solving equations.")
    residual = self.fields.get("residual")
    self._logger.stagePush("Solve")
    residual.view("RESIDUAL BEFORE SOLVE")
    self.jacobian.view()
    self.solver.solve(dispIncr, self.jacobian, residual)
    self._logger.stagePop()

    # BEGIN TEMPORARY
    #dispIncr.view("DISPINCR SOLUTION")
    #residual.view("RESIDUAL")
    # END TEMPORARY

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    # After solving, dispTBctpdt contains the displacements at time t
    # for unconstrained DOF and displacements at time t+dt at
    # constrained DOF. We add in the displacement increments (only
    # nonzero at unconstrained DOF) so that after poststep(),
    # dispTBctpdt contains the displacement field at time t+dt.
    dispIncr = self.fields.get("dispIncr(t)")
    dispTBctpdt = self.fields.solution()
    dispTBctpdt += dispIncr

    Formulation.poststep(self, t, dt)

    self._stepCount += 1
    self._logger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)

    import journal
    self._debug = journal.debug(self.name)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Implicit.
  """
  return Implicit()


# End of file 
