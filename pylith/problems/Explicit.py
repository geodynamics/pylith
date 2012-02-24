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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
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
from problems import Explicit as ModuleExplicit
from pylith.utils.profiling import resourceUsageString

# Explicit class
class Explicit(Formulation, ModuleExplicit):
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

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Formulation.Inventory):
    """
    Python object for managing Formulation facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Formulation facilities and properties.
    ##
    ## \b Properties
    ## @li \b norm_viscosity Normalized viscosity for numerical damping.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    normViscosity = pyre.inventory.float("norm_viscosity", default=0.1)
    normViscosity.meta['tip'] = "Normalized viscosity for numerical damping."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicit"):
    """
    Constructor.
    """
    Formulation.__init__(self, name)
    ModuleExplicit.__init__(self)
    self._loggingPrefix = "TSEx "
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityExplicit import ElasticityExplicit
    integrator = ElasticityExplicit()
    integrator.normViscosity(self.normViscosity)
    return integrator


  def initialize(self, dimension, normalizer):
    """
    Initialize problem for explicit time integration.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self._initialize(dimension, normalizer)

    from pylith.utils.petsc import MemoryLogger
    logger = MemoryLogger.singleton()
    logger.setDebug(0)
    logger.stagePush("Problem")

    # Allocate other fields, reusing layout from dispIncr
    if 0 == comm.rank:
      self._info.log("Creating other fields.")
    self.fields.add("disp(t-dt)", "displacement")
    self.fields.add("velocity(t)", "velocity")
    self.fields.add("acceleration(t)", "acceleration")
    self.fields.copyLayout("dispIncr(t->t+dt)")
    self._debug.log(resourceUsageString())

    # Setup fields and set to zero
    dispTmdt = self.fields.get("disp(t-dt)")
    dispTmdt.zero()
    dispT = self.fields.get("disp(t)")
    dispT.zero()
    residual = self.fields.get("residual")
    residual.zero()
    residual.createScatterMesh(residual.mesh())

    lengthScale = normalizer.lengthScale()
    timeScale = normalizer.timeScale()
    velocityScale = lengthScale / timeScale
    velocityT = self.fields.get("velocity(t)")
    velocityT.scale(velocityScale.value)
    velocityT.zero()

    accelerationScale = velocityScale / timeScale
    accelerationT = self.fields.get("acceleration(t)")
    accelerationT.scale(accelerationScale.value)
    accelerationT.zero()

    self._debug.log(resourceUsageString())
    logger.stagePop()

    if 0 == comm.rank:
      self._info.log("Creating Jacobian matrix.")
    self._setJacobianMatrixType()
    from pylith.topology.Jacobian import Jacobian
    self.jacobian = Jacobian(self.fields.solution(),
                             self.matrixType, self.blockMatrixOkay)
    self.jacobian.zero() # TEMPORARY, to get correct memory usage
    self._debug.log(resourceUsageString())

    logger.stagePush("Problem")
    if 0 == comm.rank:
      self._info.log("Initializing solver.")
    self.solver.initialize(self.fields, self.jacobian, self)
    self._debug.log(resourceUsageString())

    # Solve for increment in displacement field.
    for constraint in self.constraints:
      constraint.useSolnIncr(True)
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.useSolnIncr(True)

    logger.stagePop()
    logger.setDebug(0)
    self._eventLogger.eventEnd(logEvent)
    return


  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    logEvent = "%sprestep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    for constraint in self.constraints:
      constraint.setFieldIncr(t, t+dt, dispIncr)

    needNewJacobian = False
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.timeStep(dt)
      if integrator.needNewJacobian():
        needNewJacobian = True
    if needNewJacobian:
      self._reformJacobian(t, dt)

    self._eventLogger.eventEnd(logEvent)
    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    logEvent = "%sstep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self._reformResidual(t, dt)
    
    if 0 == comm.rank:
      self._info.log("Solving equations.")
    residual = self.fields.get("residual")
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    self.solver.solve(dispIncr, self.jacobian, residual)

    self._eventLogger.eventEnd(logEvent)
    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    # The velocity and acceleration at time t depends on the
    # displacement at time t+dt, we want to output BEFORE updating the
    # displacement fields so that the displacement, velocity, and
    # acceleration files are all at time t.
    if 0 == comm.rank:
      self._info.log("Writing solution fields.")
    for output in self.output.components():
      output.writeData(t, self.fields)
    self._writeData(t)

    # Update displacement field from time t to time t+dt.
    dispIncr = self.fields.get("dispIncr(t->t+dt)")
    dispT = self.fields.get("disp(t)")
    dispTmdt = self.fields.get("disp(t-dt)")

    dispTmdt.copy(dispT)
    dispT += dispIncr
    dispIncr.zero()

    # Complete post-step processing.
    Formulation.poststep(self, t, dt)

    self._eventLogger.eventEnd(logEvent)    
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)

    self.normViscosity = self.inventory.normViscosity
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Explicit.
  """
  return Explicit()


# End of file 
