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
    ## @li \b view_jacobian Flag to output Jacobian matrix when it is reformed.
    ##
    ## \b Facilities
    ## @li \b time_step Time step size manager.
    ## @li \b solver Algebraic solver.
    ## @li \b output Output manager associated with solution.
    ## @li \b jacobian_viewer Writer for Jacobian sparse matrix.

    import pyre.inventory

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

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="formulation"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="pde_formulation")
    ModuleFormulation.__init__(self)
    self.integratorsMesh = None
    self.integratorsSubMesh = None
    self.constraints = None
    self.jacobian = None
    self.fields = None
    self.solnName = None
    return


  def preinitialize(self, mesh, materials, boundaryConditions,
                    interfaceConditions, gravityField):
    """
    Create integrator for each element family.
    """
    self._setupLogging()
    logEvent = "%spreinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.timeStep.preinitialize()
    
    self.mesh = mesh
    self.integratorsMesh = []
    self.integratorsSubMesh = []
    self.constraints = []
    self.gravityField = gravityField

    self._setupMaterials(materials)
    self._setupBC(boundaryConditions)
    self._setupInterfaces(interfaceConditions)
    
    self._info.log("Pre-initializing output.")
    for output in self.output.components():
      output.preinitialize()

    self._logger.eventEnd(logEvent)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.timeStep.verifyConfiguration()

    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.verifyConfiguration()
    for constraint in self.constraints:
      constraint.verifyConfiguration()
    for output in self.output.components():
      output.verifyConfiguration(self.mesh)

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self, dimension, normalizer):
    """
    Create integrators for each element family.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.timeStep.initialize(normalizer)

    numTimeSteps = self.timeStep.numTimeSteps()
    totalTime = self.timeStep.totalTime

    from pylith.topology.SolutionFields import SolutionFields
    self.fields = SolutionFields(self.mesh)
    self._debug.log(resourceUsageString())

    if self.gravityField != None:
      self._info.log("Initializing gravity field.")
      self.gravityField.initialize()

    self._info.log("Initializing integrators.")
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.gravityField = self.gravityField
      integrator.initialize(totalTime, numTimeSteps, normalizer)
    ModuleFormulation.meshIntegrators(self, self.integratorsMesh)
    ModuleFormulation.submeshIntegrators(self, self.integratorsSubMesh)
    self._debug.log(resourceUsageString())

    self._info.log("Initializing constraints.")
    for constraint in self.constraints:
      constraint.initialize(totalTime, numTimeSteps, normalizer)
    self._debug.log(resourceUsageString())

    self._info.log("Setting up solution output.")
    for output in self.output.components():
      output.initialize(self.mesh, normalizer)
      output.writeInfo()
      output.open(totalTime, numTimeSteps)
    self._debug.log(resourceUsageString())

    self._info.log("Creating solution field.")
    solnName = self.solnField['name']
    self.fields.add(solnName, self.solnField['label'])
    self.fields.solutionName(solnName)
    solution = self.fields.solution()
    solution.vectorFieldType(solution.VECTOR)
    solution.newSection(solution.VERTICES_FIELD, dimension)
    for constraint in self.constraints:
      constraint.setConstraintSizes(solution)
    solution.allocate()
    for constraint in self.constraints:
      constraint.setConstraints(solution)
    solution.createScatter()
    self._debug.log(resourceUsageString())

    self._logger.eventEnd(logEvent)
    return


  def getStartTime(self):
    """
    Get start time for simulation.
    """
    return 0.0*second


  def getTotalTime(self):
    """
    Get total time for simulation.
    """
    return self.timeStep.totalTime


  def getTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    logEvent = "%stimestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    dt = self.timeStep.timeStep(self.integratorsMesh + self.integratorsSubMesh)

    self._logger.eventEnd(logEvent)
    return dt
  

  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    logEvent = "%sprestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._logger.eventEnd(logEvent)
    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    logEvent = "%sstep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    totalTime = self.timeStep.totalTime

    self._info.log("Writing solution fields.")
    for output in self.output.components():
      output.writeData(t+dt, self.fields)
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.poststep(t, dt, totalTime, self.fields)
    for constraint in self.constraints:
      constraint.poststep(t, dt, totalTime, self.fields)

    self._logger.eventEnd(logEvent)
    return


  def finalize(self):
    """
    Cleanup after time stepping.
    """
    logEvent = "%sfinalize" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Formulation finalize.")
    self._debug.log(resourceUsageString())
    for integrator in self.integratorsMesh + self.integratorsSubMesh:
      integrator.finalize()
    for constraint in self.constraints:
      constraint.finalize()
    for output in self.output.components():
      output.close()
    self._debug.log(resourceUsageString())

    self._logger.eventEnd(logEvent)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.timeStep = self.inventory.timeStep
    self.solver = self.inventory.solver
    self.output = self.inventory.output
    self.viewJacobian = self.inventory.viewJacobian
    self.jacobianViewer = self.inventory.jacobianViewer

    import journal
    self._debug = journal.debug(self.name)
    return


  def _setupMaterials(self, materials):
    """
    Setup materials as integrators.
    """
    from pylith.feassemble.Integrator import implementsIntegrator
    
    self._info.log("Pre-initializing materials.")
    self._debug.log(resourceUsageString())
    for material in materials.components():
      integrator = self.elasticityIntegrator()
      if not implementsIntegrator(integrator):
        raise TypeError, \
              "Could not use '%s' as an integrator for material '%s'. " \
              "Functionality missing." % (integrator.name, material.label)
      integrator.preinitialize(self.mesh, material)
      self.integratorsMesh.append(integrator)
      self._debug.log(resourceUsageString())

      self._info.log("Added elasticity integrator for material '%s'." % \
                     material.label)    
    return


  def _setupBC(self, boundaryConditions):
    """
    Setup boundary conditions as integrators or constraints.
    """
    from pylith.feassemble.Integrator import implementsIntegrator
    from pylith.feassemble.Constraint import implementsConstraint

    self._info.log("Pre-initializing boundary conditions.")
    self._debug.log(resourceUsageString())
    for bc in boundaryConditions.components():
      bc.preinitialize(self.mesh)
      foundType = False
      if implementsIntegrator(bc):
        foundType = True
        self.integratorsSubMesh.append(bc)
        self._info.log("Added boundary condition '%s' as an integrator." % \
                       bc.label)
      if implementsConstraint(bc):
        foundType = True
        self.constraints.append(bc)
        self._info.log("Added boundary condition '%s' as a constraint." % \
                       bc.label)
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

    self._info.log("Pre-initializing interior interfaces.")
    for ic in interfaceConditions.components():
      ic.preinitialize(self.mesh)
      foundType = False
      if implementsIntegrator(ic):
        foundType = True
        self.integratorsSubMesh.append(ic)
        self._info.log("Added interface condition '%s' as an integrator." % \
                       ic.label)
      if implementsConstraint(ic):
        foundType = True
        self.constraints.append(ic)
        self._info.log("Added interface condition '%s' as a constraint." % \
                       ic.label)
      if not foundType:
        raise TypeError, \
              "Could not determine whether interface condition '%s' is an " \
              "integrator or a constraint." % ic.name
    self._debug.log(resourceUsageString())    
    return
  

  def _reformJacobian(self, t, dt):
    """
    Reform Jacobian matrix for operator.
    """
    self._debug.log(resourceUsageString())
    self._info.log("Integrating Jacobian operator.")
    self._logger.stagePush("Reform Jacobian")

    self.updateSettings(self.jacobian, self.fields, t, dt)
    self.reformJacobian()

    self._logger.stagePop()

    if self.viewJacobian:
      self.jacobianViewer.write(self.jacobian, t)

    self._debug.log(resourceUsageString())
    return


  def _reformResidual(self, t, dt):
    """
    Reform residual vector for operator.
    """
    self._info.log("Integrating residual term in operator.")
    self._logger.stagePush("Reform Residual")

    self.updateSettings(self.jacobian, self.fields, t, dt)
    self.reformResidual()

    self._logger.stagePop()
    self._debug.log(resourceUsageString())
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
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    stages = ["Reform Jacobian",
              "Reform Residual",
              "Solve"]
    for stage in stages:
      logger.registerStage(stage)

    self._logger = logger
    return
  

  def _cleanup(self):
    """
    Deallocate PETSc and local data structures.
    """
    if not self.jacobian is None:
      self.jacobian.cleanup()
    if not self.fields is None:
      self.fields.cleanup()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Formulation.
  """
  return Formulation()


# End of file 
