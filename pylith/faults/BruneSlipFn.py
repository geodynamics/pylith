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

## @file pylith/faults/BruneSlipFn.py
##
## @brief Python object for slip time function that follows the
## integral of Brune's (1970) far-field time function.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn
from faults import BruneSlipFn as ModuleBruneSlipFn

# BruneSlipFn class
class BruneSlipFn(SlipTimeFn, ModuleBruneSlipFn):
  """
  Python object for slip time function that follows the integral of
  Brune's (1970) far-field time function.

  Inventory

  \b Properties
  @li None
  
  \b Facilities
  @li \b slip Spatial database of final slip
  @li \b slip_time Spatial database of slip initiation time
  @li \b rise_time Spatial database of rise time

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
  dbRiseTime.meta['tip'] = "Spatial database of rise time."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bruneslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    ModuleBruneSlipFn.__init__(self)
    self._loggingPrefix = "BrSF "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    ModuleBruneSlipFn.dbFinalSlip(self, self.inventory.dbSlip)
    ModuleBruneSlipFn.dbSlipTime(self, self.inventory.dbSlipTime)
    ModuleBruneSlipFn.dbRiseTime(self, self.inventory.dbRiseTime)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with BruneSlipFn.
  """
  return BruneSlipFn()


# End of file 
