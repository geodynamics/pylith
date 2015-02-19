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

## @file pylith/faults/LiuCosSlipFn.py
##
## @brief Sine/cosine slip time function from Liu, Archuleta, and Hartzell,
## BSSA, 2006 (doi:10.1785/0120060036) which has a rapid rise and then
## a gradual falloff with a finite duration.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn
from faults import LiuCosSlipFn as ModuleLiuCosSlipFn

# LiuCosSlipFn class
class LiuCosSlipFn(SlipTimeFn, ModuleLiuCosSlipFn):
  """
  Sine/cosine slip time function from Liu, Archuleta, and Hartzell,
  BSSA, 2006 (doi:10.1785/0120060036) which has a rapid rise and then
  a gradual falloff with a finite duration.

  Inventory

  \b Properties
  @li None
  
  \b Facilities
  @li \b slip Spatial database of final slip.
  @li \b slip_time Spatial database of slip initiation time.
  @li \b rise_time Spatial database of rise time (t95).

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
  
  from spatialdata.spatialdb.SimpleDB import SimpleDB
  
  dbSlip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB)
  dbSlip.meta['tip'] = "Spatial database of slip."
  
  dbSlipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB)
  dbSlipTime.meta['tip'] = "Spatial database of slip initiation time."
  
  dbRiseTime = pyre.inventory.facility("rise_time", family="spatial_database",
                                       factory=SimpleDB)
  dbRiseTime.meta['tip'] = "Spatial database of rise time (t95)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="liucosslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    ModuleLiuCosSlipFn.__init__(self)
    self._loggingPrefix = "LCSF "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    ModuleLiuCosSlipFn.dbFinalSlip(self, self.inventory.dbSlip)
    ModuleLiuCosSlipFn.dbSlipTime(self, self.inventory.dbSlipTime)
    ModuleLiuCosSlipFn.dbRiseTime(self, self.inventory.dbRiseTime)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with LiuCosSlipFn.
  """
  return LiuCosSlipFn()


# End of file 
