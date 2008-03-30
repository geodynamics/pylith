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

## @file pylith/faults/FaultCohesiveKin.py
##

## @brief Python object for a fault surface with kinematic
## (prescribed) slip implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive
from pylith.feassemble.Integrator import Integrator

# FaultCohesiveKin class
class FaultCohesiveKin(FaultCohesive, Integrator):
  """
  Python object for a fault surface with kinematic (prescribed) slip
  implemented with cohesive elements.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FaultCohesive.Inventory):
    """
    Python object for managing FaultCohesiveKin facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FaultCohesiveKin facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b eq_src Kinematic earthquake source information.
    ## @li \b output Output manager associated with fault data.

    import pyre.inventory

    from EqKinSrc import EqKinSrc
    eqsrc = pyre.inventory.facility("eq_src", family="eq_kinematic_src",
                                    factory=EqKinSrc)
    eqsrc.meta['tip'] = "Kinematic earthquake source information."

    from pylith.meshio.OutputFaultKin import OutputFaultKin
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputFaultKin)
    output.meta['tip'] = "Output manager associated with fault data."



  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivekin"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoKi "

    self.availableFields = \
        {'vertex': \
           {'info': ["normal-dir",
                     "final-slip",
                     "slip-time"],
            'data': ["slip", "traction-change"]},
         'cell': \
           {'info': [],
            'data': []}}
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    self._info.log("Pre-initializing fault '%s'." % self.label)
    FaultCohesive.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)
    assert(None != self.cppHandle)
    self.eqsrc.preinitialize()
    self.cppHandle.eqsrc = self.eqsrc.cppHandle

    if mesh.dimension() == 2:
      self.availableFields['vertex']['info'] += ["strike-dir"]
    elif mesh.dimension() == 3:
      self.availableFields['vertex']['info'] += ["strike-dir",
                                                 "dip-dir"]
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    FaultCohesive.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    self.eqsrc.verifyConfiguration()
    
    self._logger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps):
    """
    Initialize cohesive elements.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._info.log("Initializing fault '%s'." % self.label)

    self._logger.eventBegin(logEvent)
    self.eqsrc.initialize()
    FaultCohesive.initialize(self, totalTime, numTimeSteps)

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    FaultCohesive.poststep(self, t, dt, totalTime, fields)

    self._logger.eventEnd(logEvent)
    return


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      (field, fieldType) = self.cppHandle.vertexField(name,
                                                      self.mesh.cppHandle)
    else:
      (field, fieldType) = self.cppHandle.vertexField(name,
                                                     self.mesh.cppHandle,
                                                     fields.cppHandle)
    return (field, fieldType)


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      (field, fieldType) = self.cppHandle.cellField(name, self.mesh.cppHandle)
    else:
      (field, fieldType) = self.cppHandle.cellField(name, self.mesh.cppHandle,
                                                    fields.cppHandle)
    return (field, fieldType)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    FaultCohesive._configure(self)
    self.eqsrc = self.inventory.eqsrc
    self.output = self.inventory.output
    return


  def _createCppHandle(self):
    """
    Create handle to C++ FaultCohesiveKin.
    """
    if None == self.cppHandle:
      import pylith.faults.faults as bindings
      self.cppHandle = bindings.FaultCohesiveKin()
    return
    
  
# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveKin.
  """
  return FaultCohesiveKin()


# End of file 
