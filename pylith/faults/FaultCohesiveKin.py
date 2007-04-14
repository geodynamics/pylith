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
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivekin"):
    """
    Constructor.
    """
    FaultCohesive.__init__(self, name)
    import pylith.faults.faults as bindings
    self.cppHandle = bindings.FaultCohesiveKin()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    FaultCohesive._configure(self)
    return

  
# End of file 
