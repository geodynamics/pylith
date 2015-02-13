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

## @file pylith/topology/MeshRefiner.py
##
## @brief Python manager for refining mesh in parallel.
##
## Factory: mesh_refiner.

from pylith.utils.PetscComponent import PetscComponent

# MeshRefiner class
class MeshRefiner(PetscComponent):
  """
  Python manager for refining mesh in parallel.

  Factory: mesh_refiner
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="refiner"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="refiner")
    return


  def refine(self, mesh):
    """
    Refine mesh.
    """
    self._setupLogging()
    logEvent = "%srefine" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self._eventLogger.eventEnd(logEvent)
    return mesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    self._loggingPrefix = "Refin "
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Refinement")
    logger.initialize()
    events = ["refine"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_refiner():
  """
  Factory associated with MeshRefiner.
  """
  return MeshRefiner()


# End of file 
