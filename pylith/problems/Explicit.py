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
from pylith.utils.profiling import resourceUsageString

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
    self._loggingPrefix = "TSEx "
    self.solnField = {'name': "disp(t)",
                      'label': "displacement"}
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityExplicit import ElasticityExplicit
    return ElasticityExplicit()


  def initialize(self, dimension, normalizer):
    """
    Initialize problem for explicit time integration.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    Formulation.initialize(self, dimension, normalizer)

    self._info.log("Creating other fields and matrices.")
    self.fields.add("dispIncr(t->t+dt)", "displacement_increment")
    self.fields.add("disp(t-dt)", "displacement")
    self.fields.add("residual", "residual")
    self.fields.copyLayout("disp(t)")
    self.fields.solveSolnName("dispIncr(t->t+dt)")
    self._debug.log(resourceUsageString())

    # Set fields to zero
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    dispIncr.zero()
    dispT = self.fields.get("disp(t)")
    dispT.zero()
    dispTmdt = self.fields.get("disp(t-dt)")
    dispTmdt.zero()
    residual = self.fields.get("residual")
    residual.zero()
    # Create Petsc vectors for fields involved in solve
    residual.createVector()
    dispIncr.createVector()
    self._debug.log(resourceUsageString())

    self._info.log("Creating Jacobian matrix.")
    from pylith.topology.Jacobian import Jacobian
    self.jacobian = Jacobian(self.fields)
    self.jacobian.zero() # TEMPORARY, to get correct memory usage
    self._debug.log(resourceUsageString())

    self._info.log("Initializing solver.")
    self.solver.initialize(self.fields, self.jacobian, self)
    self._debug.log(resourceUsageString())

    # Solve for increment in displacement field.
    for constraint in self.constraints:
      constraint.useSolnIncr(True)
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.useSolnIncr(True)

    self._logger.eventEnd(logEvent)
    return


  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    logEvent = "%sprestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    dispTpdt = self.fields.get("disp(t+dt)")
    for constraint in self.constraints:
      constraint.setFieldIncr(t, t+dt, dispTpdt)

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

    self._reformResidual(t, dt)
    
    self._info.log("Solving equations.")
    residual = self.fields.get("residual")
    dispIncr = self.fields.get("dispIncr(t->t+dt")
    self.solver.solve(dispIncr, self.jacobian, residual)

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    disp = self.fields.get("disp(t)")
    dispTmdt = self.fields.get("disp(t-dt)")

    dispTmdt.copy(dispT)
    disp += dispIncr
    dispIncr.zero()

    Formulation.poststep(self, t, dt)

    self._logger.eventEnd(logEvent)    
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
