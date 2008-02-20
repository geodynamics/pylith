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
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.
    ## @li \b output Output manager associated with solution.

    import pyre.inventory

    from pylith.solver.SolverLinear import SolverLinear
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverLinear)
    solver.meta['tip'] = "Algebraic solver."

    from pylith.meshio.SingleOutput import SingleOutput
    output = pyre.inventory.facility("output", family="object_bin",
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

    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': ["displacements"]},
         'cell': \
           {'info': [],
            'data': []}}
    return


  def preinitialize(self, mesh, materials, boundaryConditions,
                    interfaceConditions):
    """
    Create integrator for each element family.
    """
    self._setupLogging()
    logEvent = "%spreinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    self.mesh = mesh
    self.integrators = []
    self.constraints = []

    self._setupMaterials(materials)
    self._setupBC(boundaryConditions)
    self._setupInterfaces(interfaceConditions)
    
    self._info.log("Pre-initializing output.")
    for output in self.output.bin:
      output.preinitialize()

    self._logger.eventEnd(logEvent)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    for integrator in self.integrators:
      integrator.verifyConfiguration()
    for constraint in self.constraints:
      constraint.verifyConfiguration()
    for output in self.output.bin:
      output.verifyConfiguration()

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self, dimension, totalTime, dt):
    """
    Create integrators for each element family.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    numTimeSteps = 1.0 + int(totalTime / dt)

    from pylith.topology.FieldsManager import FieldsManager
    self.fields = FieldsManager(self.mesh)
    self._debug.log(resourceUsageString())

    self._info.log("Initializing integrators.")
    for integrator in self.integrators:
      integrator.initialize(totalTime, numTimeSteps)
    self._debug.log(resourceUsageString())

    self._info.log("Initializing constraints.")
    for constraint in self.constraints:
      constraint.initialize(totalTime, numTimeSteps)
    self._debug.log(resourceUsageString())

    self._info.log("Setting up solution output.")
    for output in self.output.bin:
      output.initialize(self.mesh)
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


  def stableTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    logEvent = "%stimestep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("WARNING: Formulation::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second

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


  def poststep(self, t, dt, totalTime):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Writing solution fields.")
    for output in self.output.bin:
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
    for output in self.output.bin:
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
    self.solver = self.inventory.solver
    self.output = self.inventory.output

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
    for material in materials.bin:
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
    for bc in boundaryConditions.bin:
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
    for ic in interfaceConditions.bin:
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
    self._info.log("Reforming Jacobian of operator.")
    self._debug.log(resourceUsageString())
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(self.jacobian)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, t+dt, self.fields)
    petsc.mat_assemble(self.jacobian)
    self._debug.log(resourceUsageString())
    return


  def _reformResidual(self, t, dt):
    """
    Reform residual vector for operator.
    """
    self._info.log("Integrating residual term in operator.")
    residual = self.fields.getReal("residual")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(residual)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateResidual(residual, t, self.fields)

    self._info.log("Completing residual.")
    bindings.completeSection(self.mesh.cppHandle, residual)
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

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Formulation.
  """
  return Formulation()


# End of file 
