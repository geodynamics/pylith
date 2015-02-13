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

## @file pylith/faults/TimeHistorySlipFn.py
##
## @brief User-defined slip-time function with spatially variable
## amplitude and start time.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn
from faults import TimeHistorySlipFn as ModuleTimeHistorySlipFn

# TimeHistorySlipFn class
class TimeHistorySlipFn(SlipTimeFn, ModuleTimeHistorySlipFn):
  """
  User-defined slip-time function with spatially variable amplitude
  and start time.

  Inventory

  \b Properties
  @li None
  
  \b Facilities
  @li \b slip Spatial database of slip amplitude.
  @li \b slip_time Spatial database of slip initiation time.
  @li \b time_history Temporal database for slip time history function.

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
  
  from spatialdata.spatialdb.SimpleDB import SimpleDB
  
  dbSlip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB)
  dbSlip.meta['tip'] = "Spatial database of slip amplitude."
  
  dbSlipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB)
  dbSlipTime.meta['tip'] = "Spatial database of slip initiation time."
  
  from spatialdata.spatialdb.TimeHistory import TimeHistory
  dbTimeHistory = pyre.inventory.facility("time_history",
                                          family="temporal_database",
                                          factory=TimeHistory)
  dbTimeHistory.meta['tip'] = "Temporal database for slip evolution."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timehistoryslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    ModuleTimeHistorySlipFn.__init__(self)
    self._loggingPrefix = "THSF "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    ModuleTimeHistorySlipFn.dbAmplitude(self, self.inventory.dbSlip)
    ModuleTimeHistorySlipFn.dbSlipTime(self, self.inventory.dbSlipTime)
    ModuleTimeHistorySlipFn.dbTimeHistory(self, self.inventory.dbTimeHistory)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with TimeHistorySlipFn.
  """
  return TimeHistorySlipFn()


# End of file 
