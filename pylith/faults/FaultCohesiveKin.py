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
           {'info': ["normal_dir",
                     "final_slip",
                     "slip_time"],
            'data': ["slip"]},
         'cell': \
           {'info': [],
            'data': ["traction_change"]}}
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
    self.output.preinitialize(self)

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
    FaultCohesive.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    self.output.verifyConfiguration()
    return


  def initialize(self):
    """
    Initialize cohesive elements.
    """
    logEvent = "%sinit" % self._loggingPrefix

    self._info.log("Initializing fault '%s'." % self.label)

    self._logger.eventBegin(logEvent)
    self.eqsrc.initialize()
    FaultCohesive.initialize(self)
    self.output.initialize(self.quadrature.cppHandle)
    self.output.writeInfo()

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Writing fault data.")
    #self.output.writeData(t+dt)

    self._logger.eventEnd(logEvent)
    return


  def getVertexField(self, name):
    """
    Get vertex field.
    """
    return self.cppHandle.vertexField(name, self.mesh.cppHandle)


  def getCellField(self, name):
    """
    Get cell field.
    """
    return self.cppHandle.cellField(name, self.mesh.cppHandle)


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
    
  
  def _setupLogging(self):
    """
    Setup event logging.
    """
    Integrator._setupLogging(self)

    events = ["init"]
    for event in events:
      self._logger.registerEvent("%s%s" % (self._loggingPrefix, event))
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveKin.
  """
  return FaultCohesiveKin()


# End of file 
