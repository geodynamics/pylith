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

## @file pylith/faults/FaultCohesiveDynL.py
##

## @brief Python object for a fault surface with dynamic
## (friction) fault implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive
from pylith.feassemble.Integrator import Integrator
from faults import FaultCohesiveDynL as ModuleFaultCohesiveDynL

from pylith.utils.NullComponent import NullComponent

# FaultCohesiveDynL class
class FaultCohesiveDynL(FaultCohesive, Integrator, ModuleFaultCohesiveDynL):
  """
  Python object for a fault surface with kinematic (prescribed) slip
  implemented with cohesive elements.

  Inventory

  @class Inventory
  Python object for managing FaultCohesiveDynL facilities and properties.
  
  \b Properties
  @li None
  
  \b Facilities
  @li \b db_initial_tractions Spatial database for initial tractions.
  @li \b output Output manager associated with fault data.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  db = pyre.inventory.facility("db_initial_tractions", family="spatial_database",
                               factory=NullComponent)
  db.meta['tip'] = "Spatial database for initial tractions."

  from pylith.meshio.OutputFaultDyn import OutputFaultDyn
  output = pyre.inventory.facility("output", family="output_manager",
                                   factory=OutputFaultDyn)
  output.meta['tip'] = "Output manager associated with fault data."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivedynl"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoDy "

    self.availableFields = \
        {'vertex': \
           {'info': ["normal_dir"],
            'data': ["slip",
                     "traction"]},
         'cell': \
           {'info': [],
            'data': []}}
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    self._info.log("Pre-initializing fault '%s'." % self.label())
    FaultCohesive.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)

    ModuleFaultCohesiveDynL.quadrature(self, self.faultQuadrature)

    if mesh.dimension() == 2:
      self.availableFields['vertex']['info'] += ["strike_dir"]
    elif mesh.dimension() == 3:
      self.availableFields['vertex']['info'] += ["strike_dir",
                                                 "dip_dir"]

    if not isinstance(self.inventory.db, NullComponent):
      self.availableFields['vertex']['info'] += ["initial_traction"]
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    FaultCohesive.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    ModuleFaultCohesiveDynL.verifyConfiguration(self, self.mesh)

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize cohesive elements.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    self._info.log("Initializing fault '%s'." % self.label())

    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    
    FaultCohesive.initialize(self, totalTime, numTimeSteps, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    Integrator.poststep(self, t, dt, totalTime, fields)
    FaultCohesive.poststep(self, t, dt, totalTime, fields)

    self._eventLogger.eventEnd(logEvent)
    return


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      field = ModuleFaultCohesiveDynL.vertexField(self, name)
    else:
      field = ModuleFaultCohesiveDynL.vertexField(self, name, fields)
    return field


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      field = ModuleFaultCohesiveDynL.cellField(self, name)
    else:
      field = ModuleFaultCohesiveDynL.cellField(self, name, fields)
    return field


  def finalize(self):
    """
    Cleanup.
    """
    FaultCohesive.finalize(self)
    Integrator.finalize(self)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    FaultCohesive._configure(self)
    if not isinstance(self.inventory.db, NullComponent):
      ModuleFaultCohesiveDynL.dbInitialTract(self, self.inventory.db)
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to C++ FaultCohesiveDynL.
    """
    ModuleFaultCohesiveDynL.__init__(self)
    return
    
  
  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFault("Fault", self)
    self.perfLogger.logFields("Fault", self.fields())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveDynL.
  """
  return FaultCohesiveDynL()


# End of file 
