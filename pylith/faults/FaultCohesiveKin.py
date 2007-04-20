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

## @file pylith/faults/FaultCohesiveKin.py
##

## @brief Python object for a fault surface with kinematic
## (prescribed) slip implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive

# FaultCohesiveKin class
class FaultCohesiveKin(FaultCohesive):
  """
  Python object for a fault surface with kinematic (prescribed) slip
  implemented with cohesive elements.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FaultCohesive.Inventory):
    """
    Python object for managing FaultCohesiveKin facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FaultCohesiveKin facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b slip Spatial database of final slip
    ## @li \b slip_rate Spatial database of peak slip rate
    ## @li \b slip_time Spatial database of slip initiation time
    ## @li \b slip_function Analytical form for slip time function

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB

    slip = pyre.inventory.facility("slip", family="spatial_database",
                                   factory=SimpleDB, args=["slip"])
    slip.meta['tip'] = "Spatial database of final slip."

    slipRate = pyre.inventory.facility("slip_rate", family="spatial_database",
                                       factory=SimpleDB,
                                       args=["slip rate"])
    slipRate.meta['tip'] = "Spatial database of peak slip rate."

    slipTime = pyre.inventory.facility("slip_time", family="spatial_database",
                                       factory=SimpleDB,
                                       args=["slip time"])
    slipTime.meta['tip'] = "Spatial database of slip initiation time."

    from BruneSlipFn import BruneSlipFn
    slipFn = pyre.inventory.facility("slip_function", family="slip_time_fn",
                                     factory=BruneSlipFn)
    slipFn.meta['tip'] = "Analytical form for slip time function."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivekin"):
    """
    Constructor.
    """
    FaultCohesive.__init__(self, name)
    import pylith.faults.faults as bindings
    self.cppHandle = bindings.FaultCohesiveKin()
    return


  def initialize(self, mesh):
    """
    Initialize cohesive elements.
    """
    FaultCohesive.initialize(self, mesh)
    self.slip.initialize()
    self.slipRate.initialize()
    self.slipTime.initialize()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    FaultCohesive._configure(self)
    slip = self.inventory.slip
    slipRate = self.inventory.slipRate
    slipTime = self.inventory.slipTime
    slipFn = self.inventory.slipFn
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveKin.
  """
  return FaultCohesiveKin()


# End of file 
