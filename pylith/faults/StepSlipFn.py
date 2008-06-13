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

## @file pylith/faults/StepSlipFn.py
##
## @brief Python object for a step-function slip time function.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn

# StepSlipFn class
class StepSlipFn(SlipTimeFn):
  """
  Python object for a step-function slip time function.

  Factory: slip_time_fn
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SlipTimeFn.Inventory):
    """
    Python object for managing StepSlipFn facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing StepSlipFn facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b slip Spatial database of final slip.
    ## @li \b slip_time Spatial database of slip initiation time.

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB

    slipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB)
    slipTime.meta['tip'] = "Spatial database of slip initiation time."

    slip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB)
    slip.meta['tip'] = "Spatial database of final slip."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="stepslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    self._loggingPrefix = "StSF "
    return


  def initialize(self):
    """
    Initialize.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.slip.initialize()
    self.slipTime.initialize()
    assert(None != self.cppHandle)

    self.cppHandle.dbFinalSlip = self.slip.cppHandle
    self.cppHandle.dbSlipTime = self.slipTime.cppHandle

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
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import pylith.faults.faults as bindings
      self.cppHandle = bindings.StepSlipFn()
    return
  
  
# FACTORIES ////////////////////////////////////////////////////////////

def slip_time_fn():
  """
  Factory associated with StepSlipFn.
  """
  return StepSlipFn()


# End of file 
