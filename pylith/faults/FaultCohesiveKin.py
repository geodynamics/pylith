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

    import pyre.inventory

    from EqKinSrc import EqKinSrc
    eqsrc = pyre.inventory.facility("eq_src", family="eq_kinematic_src",
                                    factory=EqKinSrc)
    eqsrc.meta['tip'] = "Kinematic earthquake source information."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivekin"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoKi "
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
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    FaultCohesive.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
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
    self._logger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    FaultCohesive._configure(self)
    self.eqsrc = self.inventory.eqsrc
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
