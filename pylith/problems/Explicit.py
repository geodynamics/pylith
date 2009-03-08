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
    self._loggingPrefix = "TSEx "
    self.solnField = {'name': "dispT",
                      'label': "displacements"}
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
    self.fields.addReal("dispTpdt")
    self.fields.addReal("dispTmdt")
    self.fields.addReal("residual")
    self.fields.createHistory(["dispTpdt", "dispT", "dispTmdt"])    
    self.fields.copyLayout("dispT")
    self.jacobian = self.mesh.createMatrix(self.fields.getSolution())

    self.solver.initialize(self.mesh, self.fields.getSolution())

    # Solve for total displacement field
    for constraint in self.constraints:
      constraint.useSolnIncr(False)
    for integrator in self.integrators:
      integrator.useSolnIncr(False)

    self._logger.eventEnd(logEvent)
    return


  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    logEvent = "%sprestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    dispTpdt = self.fields.getReal("dispTpdt")
    for constraint in self.constraints:
      constraint.setField(t+dt, dispTpdt)

    needNewJacobian = False
    for integrator in self.integrators:
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
    residual = self.fields.getReal("residual")
    self.solver.solve(self.fields.getReal("dispTpdt"), self.jacobian, residual)

    # BEGIN TEMPORARY
    #import pylith.topology.topology as bindings
    #bindings.sectionView(residual, "RHS");
    #bindings.sectionView(self.fields.getReal("dispTpdt"), "SOLUTION");
    #import pylith.utils.petsc as petscbindings
    #print "JACOBIAN"
    #petscbindings.mat_view(self.jacobian)
    # END TEMPORARY
    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    self.fields.shiftHistory()
    if not self.solver.guessZero:
      import pylith.topology.topology as bindings
      bindings.copyRealSection(self.fields.getReal("dispTpdt"),
                               self.fields.getReal("dispT"))

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
