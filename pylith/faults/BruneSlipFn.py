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

## @file pylith/faults/BruneSlipFn.py
##
## @brief Python object for slip time function that follows the
## integral of Brune's (1970) far-field time function.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn

# BruneSlipFn class
class BruneSlipFn(SlipTimeFn):
  """
  Python object for slip time function that follows the integral of
  Brune's (1970) far-field time function.

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SlipTimeFn.Inventory):
    """
    Python object for managing BruneSlipFn facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BruneSlipFn facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b slip Spatial database of final slip
    ## @li \b slip_time Spatial database of slip initiation time
    ## @li \b slip_rate Spatial database of peak slip rate

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB

    slip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB, args=["slip"])
    slip.meta['tip'] = "Spatial database of slip."

    slipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB,
                                       args=["slip time"])
    slipTime.meta['tip'] = "Spatial database of slip initiation time."

    slipRate = pyre.inventory.facility("slip_rate", family="spatial_database",
                                       factory=SimpleDB,
                                       args=["slip rate"])
    slipRate.meta['tip'] = "Spatial database of peak slip rate."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bruneslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    import pylith.faults.faults as bindings
    self.cppHandle = bindings.BruneSlipFn()
    return


  def initialize(self):
    """
    Initialize.
    """
    SlipTimeFn.initialize(self)

    self.slip.initialize()
    self.slipTime.initialize()
    self.slipRate.initialize()

    self.cppHandle.dbFinalSlip = self.slip.cppHandle
    self.cppHandle.dbSlipTime = self.slipTime.cppHandle
    self.cppHandle.dbPeakRate = self.slipRate.cppHandle
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    self.slip = self.inventory.slip
    self.slipTime = self.inventory.slipTime
    self.slipRate = self.inventory.slipRate
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with BruneSlipFn.
  """
  return slip_time_fn()


# End of file 
