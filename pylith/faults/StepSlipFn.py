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

## @file pylith/faults/StepSlipFn.py
##
## @brief Python object for a step-function slip time function.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn
from faults import StepSlipFn as ModuleStepSlipFn

# StepSlipFn class
class StepSlipFn(SlipTimeFn, ModuleStepSlipFn):
  """
  Python object for a step-function slip time function.

  Inventory

  \b Properties
  @li None

  \b Facilities
  @li \b slip Spatial database of final slip.
  @li \b slip_time Spatial database of slip initiation time.

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  from spatialdata.spatialdb.SimpleDB import SimpleDB

  dbSlipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                     factory=SimpleDB)
  dbSlipTime.meta['tip'] = "Spatial database of slip initiation time."
  
  dbSlip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB)
  dbSlip.meta['tip'] = "Spatial database of final slip."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="stepslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    ModuleStepSlipFn.__init__(self)
    self._loggingPrefix = "StSF "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    ModuleStepSlipFn.dbSlipTime(self, self.inventory.dbSlipTime)
    ModuleStepSlipFn.dbFinalSlip(self, self.inventory.dbSlip)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with StepSlipFn.
  """
  return StepSlipFn()


# End of file 
