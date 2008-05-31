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

## @file pylith/faults/ConstRateSlipFn.py
##
## @brief Python object for a constant slip rate slip time function.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn

# ConstRateSlipFn class
class ConstRateSlipFn(SlipTimeFn):
  """
  Python object for a constant slip rate slip time function.

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SlipTimeFn.Inventory):
    """
    Python object for managing ConstRateSlipFn facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing ConstRateSlipFn facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b slip_rate Spatial database of peak slip rate
    ## @li \b slip_time Spatial database of slip initiation time

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB

    slipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB)
    slipTime.meta['tip'] = "Spatial database of slip initiation time."

    slipRate = pyre.inventory.facility("slip_rate", family="spatial_database",
                                       factory=SimpleDB)
    slipRate.meta['tip'] = "Spatial database of peak slip rate."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bruneslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    self._loggingPrefix = "CrSF "
    return


  def initialize(self):
    """
    Initialize.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.slipRate.initialize()
    self.slipTime.initialize()
    assert(None != self.cppHandle)

    self.cppHandle.dbSlipRate = self.slipRate.cppHandle
    self.cppHandle.dbSlipTime = self.slipTime.cppHandle

    self._logger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    self.slipRate = self.inventory.slipRate
    self.slipTime = self.inventory.slipTime
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import pylith.faults.faults as bindings
      self.cppHandle = bindings.ConstRateSlipFn()
    return
  
  
# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with ConstRateSlipFn.
  """
  return slip_time_fn()


# End of file 
