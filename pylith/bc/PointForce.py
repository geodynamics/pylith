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

## @file pylith/bc/PointForce.py
##
## @brief Python object for managing a point force boundary condition
## with a set of points.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from TimeDependentPoints import TimeDependentPoints
from pylith.feassemble.Integrator import Integrator
from bc import PointForce as ModulePointForce

# PointForce class
class PointForce(BoundaryCondition, 
                 TimeDependentPoints, 
                 Integrator, 
                 ModulePointForce):
  """
  Python object for managing a point force boundary condition.

  Factory: boundary_condition
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="pointforce"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "PFBC "
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self, self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize PointForce boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Initializing point forces '%s'." % self.label())

    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    BoundaryCondition.initialize(self, totalTime, numTimeSteps, normalizer)

    self._eventLogger.eventEnd(logEvent)    
    return
  

  def finalize(self):
    """
    Cleanup.
    """
    BoundaryCondition.finalize(self)
    Integrator.finalize(self)
    self._modelMemoryUse()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    TimeDependentPoints._configure(self)
    return


  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModulePointForce.__init__(self)
    return
  

  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFields("BoundaryConditions", self.parameterFields())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with PointForce.
  """
  return PointForce()

  
# End of file 
