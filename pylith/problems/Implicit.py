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
    Formulation.initialize(self, dimension, normalizer)

    # Allocate other fields, reusing layout from dispIncr
    self._info.log("Creating other fields.")
    self._info.log("Creating solution field.")
    from pylith.utils.petsc import MemoryLogger
    logger = MemoryLogger.singleton()
    #logger.setDebug(1)
    logger.stagePush("Problem")
    self.fields.copyLayout("dispIncr(t->t+dt)")

    # Setup fields and set to zero
    dispT = self.fields.get("disp(t)")
    dispT.zero()
    residual = self.fields.get("residual")
    residual.zero()
    residual.createVector()
    self._debug.log(resourceUsageString())

    self._info.log("Creating Jacobian matrix.")
    from pylith.topology.Jacobian import Jacobian
    self.jacobian = Jacobian(self.fields, self.matrixType)
    self.jacobian.zero() # TEMPORARY, to get correct memory usage
    self._debug.log(resourceUsageString())

    self._info.log("Initializing solver.")
    self.solver.initialize(self.fields, self.jacobian, self)
    self._debug.log(resourceUsageString())

    # Initial time step solves for total displacement field, not increment
    self._stepCount = 0
    for constraint in self.constraints:
      constraint.useSolnIncr(False)
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.useSolnIncr(False)

    logger.stagePop()
    logger.setDebug(0)
    return


  def getStartTime(self):
    """
    Get time at which time stepping should start.
    """
    dt = self.timeStep.timeStep(self.mesh,
                                self.integratorsMesh + self.integratorsSubMesh)
    return -dt


  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    
    # If finishing first time step, then switch from solving for total
    # displacements to solving for incremental displacements
    needNewJacobian = False
    if 1 == self._stepCount:
      self._info.log("Switching from total field solution to incremental " \
                     "field solution.")
      for constraint in self.constraints:
        constraint.useSolnIncr(True)
      for integrator in self.integratorsMesh + self.integratorsSubMesh:
        integrator.useSolnIncr(True)
      needNewJacobian = True

    self._info.log("Setting constraints.")
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    if 0 == self._stepCount:
      for constraint in self.constraints:
        constraint.setField(t+dt, dispIncr)
    else:
      for constraint in self.constraints:
        constraint.setFieldIncr(t, t+dt, dispIncr)

    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.timeStep(dt)
      if integrator.needNewJacobian():
        needNewJacobian = True
    if needNewJacobian:
      self._reformJacobian(t, dt)

    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    dispIncr.zero()

    self._reformResidual(t+dt, dt)

    self._info.log("Solving equations.")
    residual = self.fields.get("residual")
    self._logger.stagePush("Solve")
    self.solver.solve(dispIncr, self.jacobian, residual)
    self._logger.stagePop()

    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    # Update displacement field from time t to time t+dt.
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    disp = self.fields.get("disp(t)")
    disp += dispIncr
    dispIncr.zero()

    Formulation.poststep(self, t, dt)

    self._stepCount += 1
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
