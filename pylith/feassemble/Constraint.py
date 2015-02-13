#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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
              "setField",
              "poststep",
              "writeData",
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


  def poststep(self, t, dt, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    return
  

  def writeData(self, t, fields):
    """
    Write data at time t.
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
              "setSizes",
              "constraints",
              "setField",
              "poststep",
              "write",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# End of file 
