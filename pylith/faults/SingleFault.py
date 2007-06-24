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

## @file pylith/faults/SingleFault.py
##
## @brief Python faults container with one fault.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# SingleFault class
class SingleFault(ObjectBin):
  """
  Python faults container with one material.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
    """
    Python object for managing SingleFault facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing SingleFault facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b fault Fault in problem

    import pyre.inventory

    from FaultCohesiveKin import FaultCohesiveKin
    fault = pyre.inventory.facility("fault", family="fault",
                                    factory=FaultCohesiveKin)
    fault.meta['tip'] = "Fault in problem."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fault"):
    """
    Constructor.
    """
    ObjectBin.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set attributes from inventory.
    """
    ObjectBin._configure(self)
    self.bin = [self.inventory.fault]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with SingleFault.
  """
  return SingleFault()


# End of file 
