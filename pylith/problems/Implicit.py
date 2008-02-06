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
    self.solnField = {'name': "dispTBctpdt",
                      'label': "displacements"}
    self._step0 = None
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityImplicit import ElasticityImplicit
    return ElasticityImplicit()


  def initialize(self, dimension, totalTime, dt):
    """
    Initialize problem for implicit time integration.
    """
    Formulation.initialize(self, dimension, totalTime, dt)

    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    self._info.log("Creating other fields.")
    self._debug.log(resourceUsageString())
    self.fields.addReal("dispIncr")
    self.fields.addReal("residual")
    self.fields.copyLayout("dispTBctpdt")
    self._debug.log(resourceUsageString())

    self._info.log("Creating Jacobian matrix.")
    self.jacobian = self.mesh.createMatrix(self.fields.getSolution())
    self._debug.log(resourceUsageString())

    self._info.log("Initializing solver.")
    self.solver.initialize(self.mesh, self.fields.getSolution())
    self._debug.log(resourceUsageString())

    # Initial time step solves for total displacement field, not increment
    self._step0 = True
    for constraint in self.constraints:
      constraint.useSolnIncr(False)
    for integrator in self.integrators:
      integrator.useSolnIncr(False)

    self._logger.eventEnd(logEvent)
    return


  def startTime(self, dt):
    """
    Get time at which time stepping should start.
    """
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
    dispTBctpdt = self.fields.getReal("dispTBctpdt")
    for constraint in self.constraints:
      constraint.setField(t+dt, dispTBctpdt)

    needNewJacobian = False
    for integrator in self.integrators:
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

    dispIncr = self.fields.getReal("dispIncr")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(dispIncr)

    self._reformResidual(t+dt, dt)

    self._info.log("Solving equations.")
    residual = self.fields.getReal("residual")
    self.solver.solve(dispIncr, self.jacobian, residual)

    # BEGIN TEMPORARY
    import pylith.topology.topology as bindings
    bindings.sectionView(self.fields.getReal("dispIncr"), "SOLUTION");
    bindings.sectionView(self.fields.getReal("residual"), "RESIDUAL");
    # END TEMPORARY

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime):
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
    import pylith.topology.topology as bindings
    dispIncr = self.fields.getReal("dispIncr")
    disp = self.fields.getSolution()
    bindings.addRealSections(disp, disp, dispIncr)

    self._info.log("Updating integrators states.")
    for integrator in self.integrators:
      integrator.updateState(t+dt, self.fields)

    # If finishing first time step, then switch from solving for total
    # displacements to solving for incremental displacements
    if self._step0 and (t + dt) < totalTime:
      self._info.log("Switching from total field solution to incremental " \
                     "field solution.")
      for constraint in self.constraints:
        constraint.useSolnIncr(True)
      for integrator in self.integrators:
        integrator.useSolnIncr(True)
      self._reformJacobian(t, dt)
      self._step0 = False

    self._logger.eventEnd(logEvent)
    
    Formulation.poststep(self, t, dt, totalTime)
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
