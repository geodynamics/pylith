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

## @file pylith/faults/SlipTimeFn.py
##

## @brief Python abstract base class for kinematic slip time function.
##
## Factory: slip_time_fn

from pylith.utils.PetscComponent import PetscComponent

# SlipTimeFn class
class SlipTimeFn(PetscComponent):
  """
  Python abstract base class for kinematic slip time function.

  Factory: slip_time_fn
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="sliptimefn"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="sliptimefn")
    return


  def preinitialize(self):
    """
    Do pre-initialization setup.
    """
    self._setupLogging()
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self):
    """
    Initialize.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self._eventLogger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscComponent._configure(self)
    return

  
  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("Slip Time Function")
    logger.initialize()

    events = ["verify",
              "init"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# End of file 
