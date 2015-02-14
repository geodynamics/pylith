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


## @file pylith/problems/Formulation.py
##
## @brief Python abstract base class for formulations of solving equations.
##
## Factory: pde_formulation

from pylith.utils.PetscComponent import PetscComponent
from problems import Formulation as ModuleFormulation

from pylith.utils.profiling import resourceUsageString
from pyre.units.time import second

# ITEM FACTORIES ///////////////////////////////////////////////////////

def outputFactory(name):
  """
  Factory for material items.
  """
  from pyre.inventory import facility
  from pylith.meshio.OutputSoln import OutputSoln
  return facility(name, family="output_manager", factory=OutputSoln)


# Formulation class
class Formulation(PetscComponent, ModuleFormulation):
  """
  Python abstract base class for formulations of solving equations.

  In general, we use some explicit or implicit formulation of the PDEs
  to create a linear form, [A]{u}={b} that we can solve.

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Formulation facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Formulation facilities and properties.
    ##
    ## \b Properties
    ## @li \b matrix_type Type of PETSc sparse matrix.
    ## @li \b split_fields Split solution fields into displacements and Lagrange constraints.
    ## @li \b use_custom_constraint_pc Use custom preconditioner for Lagrange constraints.
    ## @li \b view_jacobian Flag to output Jacobian matrix when it is reformed.
    ##
    ## \b Facilities
    ## @li \b time_step Time step size manager.
    ## @li \b solver Algebraic solver.
    ## @li \b output Output manager associated with solution.
    ## @li \b jacobian_viewer Writer for Jacobian sparse matrix.

    import pyre.inventory

    matrixType = pyre.inventory.str("matrix_type", default="unknown")
    matrixType.meta['tip'] = "Type of PETSc sparse matrix."

    useSplitFields = pyre.inventory.bool("split_fields", default=False)
    useSplitFields.meta['tip'] = "Split solution fields into displacements "\
        "and Lagrange multipliers for separate preconditioning."

    useCustomConstraintPC = pyre.inventory.bool("use_custom_constraint_pc",
                                                default=False)
    useCustomConstraintPC.meta['tip'] = "Use custom preconditioner for " \
                                        "Lagrange constraints."

    viewJacobian = pyre.inventory.bool("view_jacobian", default=False)
    viewJacobian.meta['tip'] = "Write Jacobian matrix to binary file."
    
    from TimeStepUniform import TimeStepUniform
    timeStep = pyre.inventory.facility("time_step", family="time_step",
                                       factory=TimeStepUniform)
    timeStep.meta['tip'] = "Time step size manager."

    from SolverLinear import SolverLinear
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverLinear)
    solver.meta['tip'] = "Algebraic solver."

    from pylith.meshio.SingleOutput import SingleOutput
    output = pyre.inventory.facilityArray("output",
                                          itemFactory=outputFactory,
                                          factory=SingleOutput)
    output.meta['tip'] = "Output managers associated with solution."

    from pylith.topology.JacobianViewer import JacobianViewer
    jacobianViewer = pyre.inventory.facility("jacobian_viewer", 
                                             family="jacobian_viewer",
                                             factory=JacobianViewer)
    jacobianViewer.meta['tip'] = "Writer for Jacobian sparse matrix."

    from pylith.perf.MemoryLogger import MemoryLogger
    perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                         factory=MemoryLogger)
    perfLogger.meta['tip'] = "Performance and memory logging."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="formulation"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="formulation")
    # ModuleFormulation constructor called in base clase
    self.integrators = None
    self.constraints = None
    self.jacobian = None
    self.fields = None
    return


  def preinitialize(self, mesh, materials, boundaryConditions,
                    interfaceConditions, gravityField):
    """
    Create integrator for each element family.
    """
    self._setupLogging()
    logEvent = "%spreinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self.timeStep.preinitialize()

    import weakref
    self.mesh = weakref.ref(mesh)
    self.integrators = []
    self.constraints = []
    self.gravityField = gravityField

    self.solver.preinitialize()
    self._setupMaterials(materials)
    self._setupBC(boundaryConditions)
    self._setupInterfaces(interfaceConditions)

    if 0 == comm.rank:
      self._info.log("Pre-initializing output.")
    for output in self.output.components():
      output.preinitialize()

    self._eventLogger.eventEnd(logEvent)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self.timeStep.verifyConfiguration()

    for integrator in self.integrators:
      integrator.verifyConfiguration()
    for constraint in self.constraints:
      constraint.verifyConfiguration()
    for output in self.output.components():
      output.verifyConfiguration(self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self, dimension, normalizer):
    """
    Initialize formulation.
    """
    raise NotImplementedError("Please implement 'initialize' in derived class.")
  

  def getStartTime(self):
    """
    Get start time for simulation.
    """
    return self.timeStep.startTimeN


  def getTotalTime(self):
    """
    Get total time for simulation.
    """
    return self.timeStep.totalTimeN # Nondimensionalized total time


  def getTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    logEvent = "%stimestep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    dt = self.timeStep.timeStep(self.mesh(), self.integrators)

    self._eventLogger.eventEnd(logEvent)
    return dt
  

  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    logEvent = "%sprestep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self._eventLogger.eventEnd(logEvent)
    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    logEvent = "%sstep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self._eventLogger.eventEnd(logEvent)
    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    for integrator in self.integrators:
      integrator.poststep(t, dt, self.fields)
    for constraint in self.constraints:
      constraint.poststep(t, dt, self.fields)

    self._eventLogger.eventEnd(logEvent)
    return


  def finalize(self):
    """
    Cleanup after time stepping.
    """
    logEvent = "%sfinalize" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Formulation finalize.")
    self._debug.log(resourceUsageString())
    for integrator in self.integrators:
      integrator.finalize()
    for constraint in self.constraints:
      constraint.finalize()
    for output in self.output.components():
      output.close()
      output.finalize()
    self._debug.log(resourceUsageString())
    
    self._modelMemoryUse()

    self._eventLogger.eventEnd(logEvent)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.matrixType = self.inventory.matrixType
    self.timeStep = self.inventory.timeStep
    self.solver = self.inventory.solver
    self.output = self.inventory.output
    self.viewJacobian = self.inventory.viewJacobian
    self.jacobianViewer = self.inventory.jacobianViewer
    self.perfLogger = self.inventory.perfLogger

    import journal
    self._debug = journal.debug(self.name)

    if self.inventory.useCustomConstraintPC and \
           not self.inventory.useSplitFields:
      print "WARNING: Request to use custom preconditioner for Lagrange " \
            "constraints without splitting fields. " \
            "Setting split fields flag to 'True'."
      self.inventory.useSplitFields = True

    ModuleFormulation.splitFields(self, self.inventory.useSplitFields)
    ModuleFormulation.useCustomConstraintPC(self, self.inventory.useCustomConstraintPC)

    return


  def _setJacobianMatrixType(self):
    """
    Determine appropriate PETSc matrix type for Jacobian matrix.
    """
    # Mapping from symmetric matrix type to nonsymmetric matrix type
    matrixMap = {'sbaij': 'baij',
                 'seqsbaij': 'seqbaij',
                 'mpisbaij': 'mpibaij',
                 'unknown': 'aij'}
    isJacobianSymmetric = True
    for integrator in self.integrators:
      if not integrator.isJacobianSymmetric():
        isJacobianSymmetric = False
    if not isJacobianSymmetric:
      if self.matrixType in matrixMap.keys():
        print "WARNING: Jacobian matrix will not be symmetric.\n" \
              "         Switching matrix type from '%s' to '%s'." % \
              (self.matrixType, matrixMap[self.matrixType])
        self.matrixType = matrixMap[self.matrixType]
    self.blockMatrixOkay = True
    if self.matrixType == "unknown" and self.solver.useCUDA:
      self.matrixType = "aijcusp"
    for constraint in self.constraints:
      numDimConstrained = constraint.numDimConstrained()
      if numDimConstrained > 0 and self.mesh().dimension() != numDimConstrained:
        self.blockMatrixOkay = False
    return


  def _setupMaterials(self, materials):
    """
    Setup materials as integrators.
    """
    from pylith.feassemble.Integrator import implementsIntegrator
    
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Pre-initializing materials.")
    self._debug.log(resourceUsageString())
    for material in materials.components():
      integrator = self.elasticityIntegrator()
      if not implementsIntegrator(integrator):
        raise TypeError, \
              "Could not use '%s' as an integrator for material '%s'. " \
              "Functionality missing." % (integrator.name, material.label())
      integrator.preinitialize(self.mesh(), material)
      self.integrators.append(integrator)
      self._debug.log(resourceUsageString())

      if 0 == comm.rank:
        self._info.log("Added elasticity integrator for material '%s'." % material.label())
    return


  def _setupBC(self, boundaryConditions):
    """
    Setup boundary conditions as integrators or constraints.
    """
    from pylith.feassemble.Integrator import implementsIntegrator
    from pylith.feassemble.Constraint import implementsConstraint
    from pylith.bc.PointForce import PointForce

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Pre-initializing boundary conditions.")
    self._debug.log(resourceUsageString())
    for bc in boundaryConditions.components():
      bc.preinitialize(self.mesh())
      foundType = False
      if implementsIntegrator(bc):
        foundType = True
        self.integrators.append(bc)
        if 0 == comm.rank:
          self._info.log("Added boundary condition '%s' as an integrator." % \
                           bc.label())
      if implementsConstraint(bc):
        foundType = True
        self.constraints.append(bc)
        if 0 == comm.rank:
          self._info.log("Added boundary condition '%s' as a constraint." % \
                           bc.label())
      if not foundType:
        raise TypeError, \
              "Could not determine whether boundary condition '%s' is an " \
              "integrator or a constraint." % bc.name
    self._debug.log(resourceUsageString())    
    return


  def _setupInterfaces(self, interfaceConditions):
    """
    Setup interfaces as integrators or constraints.
    """
    from pylith.feassemble.Integrator import implementsIntegrator
    from pylith.feassemble.Constraint import implementsConstraint

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Pre-initializing interior interfaces.")
    for ic in interfaceConditions.components():
      ic.preinitialize(self.mesh())
      foundType = False
      if implementsIntegrator(ic):
        foundType = True
        self.integrators.append(ic)
        if 0 == comm.rank:
          self._info.log("Added interface condition '%s' as an integrator." % \
                           ic.label())
      if implementsConstraint(ic):
        foundType = True
        self.constraints.append(ic)
        if 0 == comm.rank:
          self._info.log("Added interface condition '%s' as a constraint." % \
                           ic.label())
      if not foundType:
        raise TypeError, \
              "Could not determine whether interface condition '%s' is an " \
              "integrator or a constraint." % ic.name
    self._debug.log(resourceUsageString())    
    return
  

  def _initialize(self, dimension, normalizer):
    """
    Create integrators for each element family.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self.timeStep.initialize(normalizer)

    numTimeSteps = self.timeStep.numTimeSteps()
    totalTime = self.timeStep.totalTime

    from pylith.topology.SolutionFields import SolutionFields
    self.fields = SolutionFields(self.mesh())
    self._debug.log(resourceUsageString())

    if 0 == comm.rank:
      self._info.log("Initializing integrators.")
    for integrator in self.integrators:
      if not self.gravityField is None:
        integrator.gravityField(self.gravityField)
      integrator.initialize(totalTime, numTimeSteps, normalizer)
    ModuleFormulation.integrators(self, self.integrators)
    self._debug.log(resourceUsageString())

    if 0 == comm.rank:
      self._info.log("Initializing constraints.")
    for constraint in self.constraints:
      constraint.initialize(totalTime, numTimeSteps, normalizer)
    self._debug.log(resourceUsageString())

    if 0 == comm.rank:
      self._info.log("Setting up solution output.")
    for output in self.output.components():
      output.initialize(self.mesh(), normalizer)
      output.writeInfo()
      output.open(totalTime, numTimeSteps)
    self._debug.log(resourceUsageString())

    # Setup fields
    if 0 == comm.rank:
      self._info.log("Creating solution field.")
    #from pylith.utils.petsc import MemoryLogger
    #memoryLogger = MemoryLogger.singleton()
    #memoryLogger.setDebug(0)
    #memoryLogger.stagePush("Problem")
    self.fields.add("dispIncr(t->t+dt)", "displacement_increment")
    self.fields.add("disp(t)", "displacement")
    self.fields.add("residual", "residual")
    self.fields.solutionName("dispIncr(t->t+dt)")

    lengthScale = normalizer.lengthScale()
    pressureScale = normalizer.pressureScale()

    solution = self.fields.get("dispIncr(t->t+dt)")
    solution.subfieldAdd("displacement", dimension, solution.VECTOR, lengthScale.value)
    solution.subfieldAdd("lagrange_multiplier", dimension, solution.VECTOR, pressureScale.value)
    solution.subfieldsSetup()
    solution.setupSolnChart()
    solution.setupSolnDof(dimension)
    # Loop over integrators to adjust DOF layout
    for integrator in self.integrators:
      integrator.setupSolnDof(solution)
    solution.vectorFieldType(solution.VECTOR)
    solution.scale(lengthScale.value)

    for constraint in self.constraints:
      constraint.setConstraintSizes(solution)
    solution.allocate()
    solution.zeroAll()
    for constraint in self.constraints:
      constraint.setConstraints(solution)
    for integrator in self.integrators:
      integrator.checkConstraints(solution)

    #memoryLogger.stagePop()

    # This also creates a global order.
    solution.createScatter(solution.mesh())

    #memoryLogger.stagePush("Problem")
    dispT = self.fields.get("disp(t)")
    dispT.vectorFieldType(dispT.VECTOR)
    dispT.scale(lengthScale.value)

    residual = self.fields.get("residual")
    residual.vectorFieldType(residual.VECTOR)
    residual.scale(lengthScale.value)

    #memoryLogger.stagePop()
    #memoryLogger.setDebug(0)
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
    self.reformJacobian()

    self._eventLogger.stagePop()

    if self.viewJacobian:
      from pylith.mpi.Communicator import Communicator
      comm = Communicator(self.mesh().comm())
      self.jacobianViewer.view(self.jacobian, t, comm)

    self._debug.log(resourceUsageString())
    return


  def _collectNeedNewJacobian(self, flagLocal):
    """
    Aggregate needNewJacobian results across processors.
    """
    if flagLocal:
      countLocal = 1
    else:
      countLocal = 0
    import pylith.mpi.mpi as mpi
    comm = self.mesh().comm()
    countAll = mpi.allreduce_scalar_int(countLocal, mpi.mpi_max(), comm.handle)
    
    return countAll > 0


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
    self.reformResidual()

    self._eventLogger.stagePop()
    self._debug.log(resourceUsageString())
    return


  def _writeData(self, t):
    """
    Write data for time t.
    """
    logEvent = "%swrite" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    for integrator in self.integrators:
      integrator.writeData(t, self.fields)
    for constraint in self.constraints:
      constraint.writeData(t, self.fields)

    self._eventLogger.eventEnd(logEvent)
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("PDE Formulation")
    logger.initialize()

    events = ["preinit",
              "verify",
              "init",
              "timestep",
              "prestep",
              "step",
              "poststep",
              "write",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    stages = ["Reform Jacobian",
              "Reform Residual",
              "Solve"]
    for stage in stages:
      logger.registerStage(stage)

    self._eventLogger = logger
    return


  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFields('Problem', self.fields)
    self.perfLogger.logJacobian('Jacobian', 'dummy')
    self.perfLogger.logGlobalOrder('GlobalOrder', 'VectorOrder',
                                   self.fields.get('residual'))
    for integrator in self.integrators:
      self.perfLogger.logQuadrature('Quadrature', integrator.quadrature())
    return


  def _cleanup(self):
    """
    Deallocate PETSc and local data structures.
    """
    if not self.jacobian is None:
      self.jacobian.cleanup()
    if not self.fields is None:
      self.fields.cleanup()
    if not self.integrators is None:
      for integrator in self.integrators:
        integrator.cleanup()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Formulation.
  """
  return Formulation()


# End of file 
