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

## @file pylith/feassemble/Constraint.py
##
## @brief Python abstract base class for constraints on operator
## actions with finite-elements.

def implementsConstraint(obj):
  """
  Check whether object implements a constraint.
  """
  result = True
  available = dir(obj)
  required = ["preinitialize",
              "verifyConfiguration",
              "initialize",
              "setConstraintSizes",
              "setConstraints",
              "useSolnIncr",
              "setField",
              "poststep",
              "finalize"]
  for attr in required:
    if not attr in available:
      result = False
  return result



# Constraint class
class Constraint(object):
  """
  Python abstract base class for constraints on operator
  actions with finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    return


  def preinitialize(self, mesh):
    """
    Setup constraint.
    """
    self._setupLogging()
    return


  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    return
  

  def finalize(self):
    """
    Cleanup.
    """
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _setupLogging(self):
    """
    Setup event logging.
    """
    if None == self._loggingPrefix:
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Constraint")
    logger.initialize()

    events = ["verify",
              "init",
              "solnIncr",
              "setSizes",
              "constraints",
              "setField",
              "poststep",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# End of file 
