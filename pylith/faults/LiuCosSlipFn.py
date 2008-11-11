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
## @brief Sine/cosine slip time function from Liu, Archuleta, and Hartzell,
## BSSA, 2006 (doi:10.1785/0120060036) which has a rapid rise and then
## a gradual falloff with a finite duration.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn

# LiuCosSlipFn class
class LiuCosSlipFn(SlipTimeFn):
  """
  Sine/cosine slip time function from Liu, Archuleta, and Hartzell,
  BSSA, 2006 (doi:10.1785/0120060036) which has a rapid rise and then
  a gradual falloff with a finite duration.

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SlipTimeFn.Inventory):
    """
    Python object for managing LiuCosSlipFn facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing LiuCosSlipFn facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b slip Spatial database of final slip.
    ## @li \b slip_time Spatial database of slip initiation time.
    ## @li \b rise_time Spatial database of rise time (t95).

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB

    slip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB)
    slip.meta['tip'] = "Spatial database of slip."

    slipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB)
    slipTime.meta['tip'] = "Spatial database of slip initiation time."

    riseTime = pyre.inventory.facility("rise_time", family="spatial_database",
                                       factory=SimpleDB)
    riseTime.meta['tip'] = "Spatial database of rise time (t95)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bruneslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    self._loggingPrefix = "BrSF "
    return


  def initialize(self):
    """
    Initialize.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.slip.initialize()
    self.slipTime.initialize()
    self.riseTime.initialize()
    assert(None != self.cppHandle)

    self.cppHandle.dbFinalSlip = self.slip.cppHandle
    self.cppHandle.dbSlipTime = self.slipTime.cppHandle
    self.cppHandle.dbRiseTime = self.riseTime.cppHandle

    self._logger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    self.slip = self.inventory.slip
    self.slipTime = self.inventory.slipTime
    self.riseTime = self.inventory.riseTime
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import pylith.faults.faults as bindings
      self.cppHandle = bindings.LiuCosSlipFn()
    return
  
  
# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with LiuCosSlipFn.
  """
  return LiuCosSlipFn()


# End of file 
