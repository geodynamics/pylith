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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/ExplicitLumped.py
##
## @brief Python ExplicitLumped object for solving equations using an
## explicit formulation with a lumped Jacobian matrix that is stored
## as a Field.
##
## Factory: pde_formulation

from Explicit import Explicit
from problems import Explicit as ModuleExplicit
from pylith.utils.profiling import resourceUsageString

# ExplicitLumped class
class ExplicitLumped(Explicit, ModuleExplicit):
  """
  Python ExplicitLumped object for solving equations using an explicit
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

  class Inventory(Explicit.Inventory):
    """
    Python object for managing ExplicitLumped facilities and properties.

    Provide appropriate solver for lumped Jacobian as the default.
    """

    ## @class Inventory
    ## Python object for managing ExplicitLumped facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.

    import pyre.inventory

    from SolverLumped import SolverLumped
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverLumped)
    solver.meta['tip'] = "Algebraic solver."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicitlumped"):
    """
    Constructor.
    """
    Explicit.__init__(self, name)
    return


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
      self._info.log("Creating lumped Jacobian matrix.")
    from pylith.topology.topology import MeshField
    jacobian = MeshField(self.mesh)
    jacobian.newSection(jacobian.VERTICES_FIELD, dimension)
    jacobian.allocate()
    jacobian.label("jacobian")
    jacobian.vectorFieldType(jacobian.VECTOR)
    self.jacobian = jacobian
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


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Explicit._configure(self)
    self.solver = self.inventory.solver
    return


  def _reformResidual(self, t, dt):
    """
    Reform residual vector for operator.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Integrating residual term in operator.")
    self._eventLogger.stagePush("Reform Residual")

    self.updateSettings(self.jacobian, self.fields, t, dt)
    self.reformResidualLumped()

    self._eventLogger.stagePop()
    self._debug.log(resourceUsageString())
    return


  def _reformJacobian(self, t, dt):
    """
    Reform Jacobian matrix for operator.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self._debug.log(resourceUsageString())
    if 0 == comm.rank:
      self._info.log("Integrating Jacobian operator.")
    self._eventLogger.stagePush("Reform Jacobian")

    self.updateSettings(self.jacobian, self.fields, t, dt)
    self.reformJacobianLumped()

    self._eventLogger.stagePop()

    if self.viewJacobian:
      self.jacobian.view("Lumped Jacobian")

    self._debug.log(resourceUsageString())
    return


  def _cleanup(self):
    """
    Deallocate PETSc and local data structures.
    """
    if not self.fields is None:
      self.fields.cleanup()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with ExplicitLumped.
  """
  return ExplicitLumped()


# End of file 
