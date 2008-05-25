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

## @file pylith/faults/EqKinSrc.py
##

## @brief Python object for managing parameters for a kinematic
## earthquake sources.
##
## EqKinSrc is responsible for providing the value of slip at time t
## over a fault surface.
##
## Factory: eq_kinematic_src

from pyre.components.Component import Component

# EqKinSrc class
class EqKinSrc(Component):
  """
  Python object for managing parameters for a kinematic earthquake sources.

  Factory: eq_kinematic_src
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing EqKinSrc facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing EqKinSrc facilities and properties.
    ##
    ## \b Properties
    ## @li \b origin_time Origin time for earthquake rupture.
    ##
    ## \b Facilities
    ## @li \b slip_function Slip time history function.

    import pyre.inventory

    from pyre.units.time import second
    originTime = pyre.inventory.dimensional("origin_time", default=0.0*second)
    originTime.meta['tip'] = "Origin time for earthquake rupture."

    from BruneSlipFn import BruneSlipFn
    slipfn = pyre.inventory.facility("slip_function", family="slip_time_fn",
                                     factory=BruneSlipFn)
    slipfn.meta['tip'] = "Slip time history function."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="eqkinsrc"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="eqkinsrc")
    self.cppHandle = None
    self._loggingPrefix = "EqKi "
    return


  def preinitialize(self):
    """
    Do pre-initialization setup.
    """
    self._setupLogging()
    self._createCppHandle()
    self.cppHandle.originTime = self.originTime.value
    self.slipfn.preinitialize()
    self.cppHandle.slipfn = self.slipfn.cppHandle
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.slipfn.verifyConfiguration()

    self._logger.eventEnd(logEvent)
    return


  def initialize(self):
    """
    Initialize.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.slipfn.initialize()

    self._logger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.originTime = self.inventory.originTime
    self.slipfn = self.inventory.slipfn
    return

  
  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.faults.faults as bindings
      self.cppHandle = bindings.EqKinSrc()
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("FE Constraint")
    logger.initialize()

    events = ["verify",
              "init"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
  """
  Factory associated with EqKinSrc.
  """
  return EqKinSrc()


# End of file 
