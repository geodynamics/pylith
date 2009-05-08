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

## @file pylith/problems/Problem.py
##
## @brief Python abstract base class for crustal dynamics problems.
##
## Factory: problem.

from pylith.utils.PetscComponent import PetscComponent

# ITEM FACTORIES ///////////////////////////////////////////////////////

def materialFactory(name):
  """
  Factory for material items.
  """
  from pyre.inventory import facility
  from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D
  return facility(name, family="material", factory=ElasticIsotropic3D)


def bcFactory(name):
  """
  Factory for boundary condition items.
  """
  from pyre.inventory import facility
  from pylith.bc.DirichletBC import DirichletBC
  return facility(name, family="boundary_condition", factory=DirichletBC)


def faultFactory(name):
  """
  Factory for fault items.
  """
  from pyre.inventory import facility
  from pylith.faults.FaultCohesiveKin import FaultCohesiveKin
  return facility(name, family="fault", factory=FaultCohesiveKin)


# Problem class
class Problem(PetscComponent):
  """
  Python abstract base class for crustal dynamics problems.

  Factory: problem.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Problem facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Problem facilities and properties.
    ##
    ## \b Properties
    ## @li \b dimension Spatial dimension of problem space.
    ## @li \b useGravity Gravity on (true) or off (false).
    ##
    ## \b Facilities
    ## @li \b normalizer Nondimensionalizer for problem.
    ## @li \b materials Materials in problem.
    ## @li \b bc Boundary conditions.
    ## @li \b interfaces Interior surfaces with constraints or
    ##   constitutive models.
    ## @li \b gravityField Gravity field for problem (SpatialDB).

    import pyre.inventory
    from pylith.utils.EmptyBin import EmptyBin

    useGravity = pyre.inventory.bool("use_gravity", default=False)
    useGravity.meta['tip'] = "Use gravitational body forces in problem."
    
    dimension = pyre.inventory.int("dimension", default=3,
                                   validator=pyre.inventory.choice([1,2,3]))
    dimension.meta['tip'] = "Spatial dimension of problem space."

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = pyre.inventory.facility("normalizer",
                                         family="nondimensional",
                                         factory=Nondimensional)
    normalizer.meta['tip'] = "Nondimensionalizer for problem."

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facilityArray("materials",
                                             itemFactory=materialFactory,
                                             factory=Homogeneous)
    materials.meta['tip'] = "Materials in problem."

    bc = pyre.inventory.facilityArray("bc",
                                      itemFactory=bcFactory,
                                      factory=EmptyBin)
    bc.meta['tip'] = "Boundary conditions."

    interfaces = pyre.inventory.facilityArray("interfaces",
                                              itemFactory=faultFactory,
                                              factory=EmptyBin)
    interfaces.meta['tip'] = "Interior surfaces with constraints or " \
                             "constitutive models."

    from spatialdata.spatialdb.GravityField import GravityField
    gravityField = pyre.inventory.facility("gravity_field",
                                          factory=GravityField,
                                          family="spatial_database")
    gravityField.meta['tip'] = "Database used for gravity field."

    from pylith.perf.MemoryLogger import MemoryLogger
    perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                         factory=MemoryLogger)
    perfLogger.meta['tip'] = "Performance and memory logging."

  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="problem"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="problem")
    self.mesh = None
    return


  def preinitialize(self, mesh):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    raise NotImplementedError, "initialize() not implemented."
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix    
    self._logger.eventBegin(logEvent)
    
    self._info.log("Verifying compatibility of problem configuration.")
    if self.dimension != self.mesh.dimension():
      raise ValueError, \
            "Spatial dimension of problem is '%d' but mesh contains cells " \
            "for spatial dimension '%d'." % \
            (self.dimension, self.mesh.dimension())

    # Check to make sure ids of materials and interfaces are unique
    materialIds = {}
    for material in self.materials.components():
      if material.quadrature.spaceDim() != self.dimension:
        raise ValueError, \
              "Spatial dimension of problem is '%d' but quadrature " \
              "for material '%s' is for spatial dimension '%d'." % \
              (self.dimension, material.label(), material.quadrature.spaceDim())
      if material.id() in materialIds.keys():
        raise ValueError, \
            "ID values for materials '%s' and '%s' are both '%d'. " \
            "Material id values must be unique." % \
            (material.label(), materialIds[material.id()], material.id())
      materialIds[material.id()] = material.label()
      self.perfLogger.logMaterial("Problem Verification", material)
    
    for interface in self.interfaces.components():
      if interface.id() in materialIds.keys():
        raise ValueError, \
            "ID values for material '%s' and interface '%s' are both '%d'. " \
            "Material and interface id values must be unique." % \
            (materialIds[interface.id()], interface.label(), interface.id())
      materialIds[interface.id()] = interface.label()

    # Check to make sure material-id for each cell matches the id of a material
    import numpy
    idValues = numpy.array(materialIds.keys(), dtype=numpy.int32)
    self.mesh.checkMaterialIds(idValues)

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self):
    """
    Initialize integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    raise NotImplementedError, "initialize() not implemented."
    return


  def run(self, app):
    """
    Solve the problem.
    """
    raise NotImplementedError, "run() not implemented."
    return


  def finalize(self, mesh):
    """
    Cleanup.
    """
    raise NotImplementedError, "finalize() not implemented."
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    raise NotImplementedError, "checkpoint() not implemented."
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.normalizer = self.inventory.normalizer
    self.dimension = self.inventory.dimension
    self.materials = self.inventory.materials
    self.bc = self.inventory.bc
    self.interfaces = self.inventory.interfaces
    if self.inventory.useGravity:
      self.gravityField = self.inventory.gravityField
    else:
      self.gravityField = None
    self.perfLogger = self.inventory.perfLogger
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("Problem")
    logger.initialize()

    events = ["preinit",
              "verify",
              "init",
              "run",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with Problem.
  """
  return Problem()


# End of file 
