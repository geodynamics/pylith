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

## @file pylith/faults/EqKinSrc.py
##

## @brief Python object for managing parameters for a kinematic
## earthquake sources.
##
## EqKinSrc is responsible for providing the value of slip at time t
## over a fault surface.
##
## Factory: eq_kinematic_src

from pylith.utils.PetscComponent import PetscComponent
from faults import EqKinSrc as ModuleEqKinSrc

# EqKinSrc class
class EqKinSrc(PetscComponent, ModuleEqKinSrc):
  """
  Python object for managing parameters for a kinematic earthquake sources.

  Inventory

  \b Properties
  @li \b origin_time Origin time for earthquake rupture.
  
  \b Facilities
  @li \b slip_function Slip time history function.

  Factory: eq_kinematic_src
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
  
  from pyre.units.time import second
  originTime = pyre.inventory.dimensional("origin_time", default=0.0*second)
  originTime.meta['tip'] = "Origin time for earthquake rupture."
  
  from StepSlipFn import StepSlipFn
  slipfn = pyre.inventory.facility("slip_function", family="slip_time_fn",
                                   factory=StepSlipFn)
  slipfn.meta['tip'] = "Slip time history function."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="eqkinsrc"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="eqkinsrc")
    self._createModuleObj()
    self._loggingPrefix = "EqKi "
    return


  def preinitialize(self):
    """
    Do pre-initialization setup.
    """
    self._setupLogging()
    self.slipfn.preinitialize()
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self.slipfn.verifyConfiguration()

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self):
    """
    Initialize.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self.slipfn.initialize()

    self._eventLogger.eventEnd(logEvent)
    return


  def finalize(self):
    """
    Cleanup.
    """
    self.slipfn.finalize()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscComponent._configure(self)
    ModuleEqKinSrc.originTime(self, self.inventory.originTime.value)
    ModuleEqKinSrc.slipfn(self, self.inventory.slipfn)
    return

  
  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModuleEqKinSrc.__init__(self)
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("Kinematic Earthquake Source")
    logger.initialize()

    events = ["verify",
              "init"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
  """
  Factory associated with EqKinSrc.
  """
  return EqKinSrc()


# End of file 
