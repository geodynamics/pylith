#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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
from problems import Implicit as ModuleImplicit
from pylith.utils.profiling import resourceUsageString

# Implicit class
class Implicit(Formulation, ModuleImplicit):
  """
  Python Implicit object for solving equations using an implicit
  formulation.

  The formulation has the general form,

  [A(t+dt)] {u(t+dt)} = {b(t+dt)}. 

  We know the solution at time t, so we write {u(t+dt)} as {u(t)} +
  {du(t)}, where {du(t)} is the increment in the solution from time t
  to time t+dt. Thus, we solve

  [A(t+dt)] {du(t)} = {b(t+dt)} - [A(t+dt)]{u(t)}.

  We solve this system by forming the Jacobian, A, and the residual

  {r(t+dt)} = {b(t+dt)} - [A(t+dt)]{u(t)} - [A(t+dt)]{du(t)}

  which we combine into

  {r(t+dt)} = {b(t+dt)} - [A(t+dt)]{u(t)+du(t)}.

  The method reformJacobian() computes [A(t+dt)] and the method
  reformResidual computes {r(t+dt)}. Note that in forming the residual
  we compute the action [A(t+dt)]{u(t)+du(t)} and do not perform a
  matrix-vector multiplication.

  [A(t+dt)] generally depends on {u(t+dt)} as well as the current
  stresses and additional state variables.  

  For linear elastic or viscoelastic problems with constant time step
  size, A is a constant (after the elastic solution).  {b(t+dt)}
  generally depends on the loads applied for time step t+dt (including
  the contributions to the internal force vector from
  displacement/velocity BC) as well as the internal force vector
  computed from the current stresses.

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
    ModuleImplicit.__init__(self)
    self._loggingPrefix = "TSIm "
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityImplicit import ElasticityImplicit
    integrator = ElasticityImplicit()
    return integrator


  def initialize(self, dimension, normalizer):
    """
    Initialize problem for implicit time integration.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    comm = self.mesh().comm()

    self._initialize(dimension, normalizer)

    #from pylith.utils.petsc import MemoryLogger
    #memoryLogger = MemoryLogger.singleton()
    #memoryLogger.setDebug(0)
    #memoryLogger.stagePush("Problem")

    # Allocate other fields, reusing layout from dispIncr
    if 0 == comm.rank:
      self._info.log("Creating other fields.")
    self.fields.add("velocity(t)", "velocity")
    self.fields.copyLayout("dispIncr(t->t+dt)")

    # Setup fields and set to zero
    dispT = self.fields.get("disp(t)")
    dispT.zeroAll()
    residual = self.fields.get("residual")
    residual.zeroAll()
    residual.createScatter(residual.mesh())

    lengthScale = normalizer.lengthScale()
    timeScale = normalizer.timeScale()
    velocityScale = lengthScale / timeScale
    velocityT = self.fields.get("velocity(t)")
    velocityT.scale(velocityScale.value)
    velocityT.zeroAll()

    self._debug.log(resourceUsageString())
    #memoryLogger.stagePop()

    # Allocates memory for nonzero pattern and Jacobian
    if 0 == comm.rank:
      self._info.log("Creating Jacobian matrix.")
    self._setJacobianMatrixType()
    from pylith.topology.Jacobian import Jacobian
    self.jacobian = Jacobian(self.fields.solution(),
                             self.matrixType, self.blockMatrixOkay)
    self.jacobian.zero() # TEMPORARY, to get correct memory usage
    self._debug.log(resourceUsageString())

    #memoryLogger.stagePush("Problem")
    if 0 == comm.rank:
      self._info.log("Initializing solver.")
    self.solver.initialize(self.fields, self.jacobian, self)
    self._debug.log(resourceUsageString())

    #memoryLogger.stagePop()
    #memoryLogger.setDebug(0)
    return


  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    comm = self.mesh().comm()
    
    if 0 == comm.rank:
      self._info.log("Setting constraints.")
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    dispIncr.zeroAll()
    for constraint in self.constraints:
      constraint.setFieldIncr(t, t+dt, dispIncr)

    needNewJacobian = False
    for integrator in self.integrators:
      integrator.timeStep(dt)
      if integrator.needNewJacobian():
        needNewJacobian = True
    if self._collectNeedNewJacobian(needNewJacobian):
      self._reformJacobian(t, dt)

    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    comm = self.mesh().comm()

    dispIncr = self.fields.get("dispIncr(t->t+dt)")

    self._reformResidual(t+dt, dt)

    if 0 == comm.rank:
      self._info.log("Solving equations.")
    self._eventLogger.stagePush("Solve")

    residual = self.fields.get("residual")
    #self.jacobian.view() # TEMPORARY
    self.solver.solve(dispIncr, self.jacobian, residual)
    #dispIncr.view("DISP INCR") # TEMPORARY

    # DEBUGGING Verify solution makes residual 0
    #self._reformResidual(t+dt, dt)
    #residual.view("RESIDUAL")
    
    self._eventLogger.stagePop()

    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    comm = self.mesh().comm()

    # Update displacement field from time t to time t+dt.
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    disp = self.fields.get("disp(t)")
    disp.add(dispIncr)
    dispIncr.zeroAll()

    # Complete post-step processing, then write data.
    Formulation.poststep(self, t, dt)

    # Write data. Velocity at time t will be based upon displacement
    # at time t-dt and t.
    if 0 == comm.rank:
      self._info.log("Writing solution fields.")
    for output in self.output.components():
      output.writeData(t+dt, self.fields)
    self._writeData(t+dt)

    return


  def prestepElastic(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    comm = self.mesh().comm()
    
    if 0 == comm.rank:
      self._info.log("Setting constraints.")
    disp = self.fields.get("dispIncr(t->t+dt)")
    disp.zeroAll()
    for constraint in self.constraints:
      constraint.setField(t+dt, disp)

    needNewJacobian = False
    for integrator in self.integrators:
      integrator.timeStep(dt)
      if integrator.needNewJacobian():
        needNewJacobian = True
    if self._collectNeedNewJacobian(needNewJacobian):
      self._reformJacobian(t, dt)

    return


  def finalize(self):
    """
    Cleanup after time stepping.
    """
    Formulation.finalize(self)
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
