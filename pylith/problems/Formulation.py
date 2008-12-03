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

from pyre.components.Component import Component

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
class Formulation(Component):
  """
  Python abstract base class for formulations of solving equations.

  In general, we use some explicit or implicit formulation of the PDEs
  to create a linear form, [A]{u}={b} that we can solve.

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Formulation facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Formulation facilities and properties.
    ##
    ## \b Properties
    ## @li \b view_jacobian Flag to output Jacobian matrix when it is reformed.
    ## @li \b jacobian_filename Filename for Jacobian matrix.
    ## @li \b jacobian_time_width String length of time stamp in Jacobian filename.
    ## @li \b jacobian_time_constant Value used to normalize time stamp in filename.

    ##
    ## \b Facilities
    ## @li \b time_step Time step size manager.
    ## @li \b solver Algebraic solver.
    ## @li \b output Output manager associated with solution.

    import pyre.inventory

    viewJacobian = pyre.inventory.bool("view_jacobian", default=False)
    viewJacobian.meta['tip'] = "Write Jacobian matrix to binary file."
    
    jacobianFilename = pyre.inventory.str("jacobian_filename",
                                          default="jacobian.mat")
    jacobianFilename.meta['tip'] = "Filename for Jacobian matrix."

    jacobianTimeWidth = pyre.inventory.int("jacobian_time_width", default=5)
    jacobianTimeWidth.meta['tip'] = "String length of time stamp in Jacobian filename."

    jacobianTimeConstant = pyre.inventory.dimensional("jacobian_time_constant",
                                                      default=1.0*second,
                                                      validator=pyre.inventory.greater(0.0*second))
    jacobianTimeConstant.meta['tip'] = "Values used to normalize time stamp in Jacobian filename."
    
    from TimeStepUniform import TimeStepUniform
    timeStep = pyre.inventory.facility("time_step", family="time_step",
                                       factory=TimeStepUniform)
    timeStep.meta['tip'] = "Time step size manager."

    from pylith.solver.SolverLinear import SolverLinear
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverLinear)
    solver.meta['tip'] = "Algebraic solver."

    from pylith.meshio.SingleOutput import SingleOutput
    output = pyre.inventory.facilityArray("output",
                                          itemFactory=outputFactory,
                                          factory=SingleOutput)
    output.meta['tip'] = "Output managers associated with solution."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="formulation"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="pde_formulation")
    self.integrators = None
    self.constraints = None
    self.fields = None
    self.solnField = None
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
    self.integrators = []
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

    for integrator in self.integrators:
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

    from pylith.topology.FieldsManager import FieldsManager
    self.fields = FieldsManager(self.mesh)
    self._debug.log(resourceUsageString())

    if self.gravityField != None:
      self._info.log("Initializing gravity field.")
      self.gravityField.initialize()

    self._info.log("Initializing integrators.")
    for integrator in self.integrators:
      integrator.gravityField = self.gravityField
      integrator.initialize(totalTime, numTimeSteps, normalizer)
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
    self.fields.addReal(solnName)
    self.fields.solutionField(solnName)
    self.fields.setFiberDimension(solnName, dimension)
    for constraint in self.constraints:
      constraint.setConstraintSizes(self.fields.getSolution())
    self.fields.allocate(solnName)
    for constraint in self.constraints:
      constraint.setConstraints(self.fields.getSolution())
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

    dt = self.timeStep.timeStep(self.integrators)

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
    for integrator in self.integrators:
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
    for integrator in self.integrators:
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
    Component._configure(self)
    self.timeStep = self.inventory.timeStep
    self.solver = self.inventory.solver
    self.output = self.inventory.output
    self.viewJacobian = self.inventory.viewJacobian
    self.jacobianFilename = self.inventory.jacobianFilename
    self.jacobianTimeWidth = self.inventory.jacobianTimeWidth
    self.jacobianTimeConstant = self.inventory.jacobianTimeConstant

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
      self.integrators.append(integrator)
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
        self.integrators.append(bc)
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
        self.integrators.append(ic)
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
    self._logger.stagePush("Reform Jacobian")
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(self.jacobian)

    # Add in contributions that do not require assembly
    self._info.log("Reforming assembled portion of Jacobian of operator.")
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobianAssembled(self.jacobian, t+dt,
                                            self.fields)
      self._debug.log(resourceUsageString()) # TEMPORARY
    self._info.log("Flushing assembly of Jacobian of operator.")
    petsc.mat_assemble(self.jacobian, "flush_assembly")

    # Add in contributions that require assembly
    self._info.log("Reforming unassembled portion of Jacobian of operator.")
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, t+dt, self.fields)
      self._debug.log(resourceUsageString()) # TEMPORARY
    self._info.log("Doing final assembly of Jacobian of operator.")
    petsc.mat_assemble(self.jacobian, "final_assembly")
    self._logger.stagePop()

    if self.viewJacobian:
      filename = self._createJacobianFilename(t+dt)
      petsc.mat_view_binary(self.jacobian, filename)
    self._debug.log(resourceUsageString())
    return


  def _createJacobianFilename(self, t):
    """
    Create filename by extracting basename and adding a time stamp.
    """
    base = self.jacobianFilename.lstrip().rstrip()
    baseLen = len(base)
    if base.endswith(".mat"):
      basename = base[0:baseLen-4]
    else:
      basename = base
    time = int(t.value/self.jacobianTimeConstant.value)
    timeStamp = repr(time).rjust(self.jacobianTimeWidth, '0')
    filename = basename + "_t" + timeStamp + ".mat"
    return filename

      
  def _reformResidual(self, t, dt):
    """
    Reform residual vector for operator.
    """
    self._info.log("Integrating residual term in operator.")
    self._logger.stagePush("Reform Residual")
    residual = self.fields.getReal("residual")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(residual)

    # Add in contributions that require assembly
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateResidual(residual, t, self.fields)
    self._info.log("Completing residual.")
    bindings.completeSection(self.mesh.cppHandle, residual)

    # Add in contributions that do not require assembly
    self._info.log("Integrating assembled residual term in operator.")
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateResidualAssembled(residual, t, self.fields)

    self._logger.stagePop()
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("PDE Formulation")
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
  

# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Formulation.
  """
  return Formulation()


# End of file 
