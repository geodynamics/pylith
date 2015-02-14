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

## @file pylith/faults/ConstRateSlipFn.py
##
## @brief Python object for a constant slip rate slip time function.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn
from faults import ConstRateSlipFn as ModuleConstRateSlipFn

# ConstRateSlipFn class
class ConstRateSlipFn(SlipTimeFn, ModuleConstRateSlipFn):
  """
  Python object for a constant slip rate slip time function.

  Inventory

  \b Properties
  @li None
  
  \b Facilities
  @li \b slip_rate Spatial database of slip rate
  @li \b slip_time Spatial database of slip initiation time

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  from spatialdata.spatialdb.SimpleDB import SimpleDB
  
  dbSlipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB)
  dbSlipTime.meta['tip'] = "Spatial database of slip initiation time."
  
  dbSlipRate = pyre.inventory.facility("slip_rate", family="spatial_database",
                                       factory=SimpleDB)
  dbSlipRate.meta['tip'] = "Spatial database of slip rate."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="constrateslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    ModuleConstRateSlipFn.__init__(self)
    self._loggingPrefix = "CrSF "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    ModuleConstRateSlipFn.dbSlipRate(self, self.inventory.dbSlipRate)
    ModuleConstRateSlipFn.dbSlipTime(self, self.inventory.dbSlipTime)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with ConstRateSlipFn.
  """
  return ConstRateSlipFn()


# End of file 
