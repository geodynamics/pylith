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

# ITEM FACTORIES ///////////////////////////////////////////////////////

def eqsrcFactory(name):
  """
  Factory for earthquake source items.
  """
  from pyre.inventory import facility
  from EqKinSrc import EqKinSrc
  return facility(name, family="eq_kinematic_src", factory=EqKinSrc)


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
    ## @li \b eq_srcs Kinematic earthquake sources information.
    ## @li \b output Output manager associated with fault data.

    import pyre.inventory

    from SingleRupture import SingleRupture
    eqsrcs = pyre.inventory.facilityArray("eq_srcs", itemFactory=eqsrcFactory,
                                          factory=SingleRupture)
    eqsrcs.meta['tip'] = "Kinematic earthquake sources information."

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
           {'info': ["normal_dir",
                     "final_slip",
                     "slip_time"],
            'data': ["slip",
                     "traction_change"]},
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
    for eqsrc in self.eqsrcs.components():
      eqsrc.preinitialize()
    self.cppHandle.eqsrcs(self.eqsrcs.inventory.facilityNames(),
                          self.eqsrcs.components())

    for name in self.eqsrcs.inventory.facilityNames():
      self.availableFields['vertex']['info'] += ["final_slip_%s" % name]
      self.availableFields['vertex']['info'] += ["slip_time_%s" % name]

    if mesh.dimension() == 2:
      self.availableFields['vertex']['info'] += ["strike_dir"]
    elif mesh.dimension() == 3:
      self.availableFields['vertex']['info'] += ["strike_dir",
                                                 "dip_dir"]
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    FaultCohesive.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    for eqsrc in self.eqsrcs.components():
      eqsrc.verifyConfiguration()
    
    self._logger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps):
    """
    Initialize cohesive elements.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    self._info.log("Initializing fault '%s'." % self.label)

    Integrator.initialize(self, totalTime, numTimeSteps)
    
    for eqsrc in self.eqsrcs.components():
      eqsrc.initialize()
    FaultCohesive.initialize(self, totalTime, numTimeSteps)

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    Integrator.poststep(self, t, dt, totalTime, fields)
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
    self.eqsrcs = self.inventory.eqsrcs
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
