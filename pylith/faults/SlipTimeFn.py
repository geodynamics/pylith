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

## @file pylith/faults/SlipTimeFn.py
##

## @brief Python abstract base class for kinematic slip time function.
##
## Factory: slip_time_fn

from pyre.components.Component import Component

# SlipTimeFn class
class SlipTimeFn(Component):
  """
  Python abstract base class for kinematic slip time function.

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

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB

    slip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB, args=["slip"])
    slip.meta['tip'] = "Spatial database of slip."

    slipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB,
                                       args=["slip time"])
    slipTime.meta['tip'] = "Spatial database of slip initiation time."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="sliptimefn"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="sliptimefn")
    self.cppHandle = None
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.slip = self.inventory.slip
    self.slipTime = self.inventory.slipTime
    return

  
# End of file 
