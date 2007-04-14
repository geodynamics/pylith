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
## Factory: faults_bin

from FaultsBin import FaultsBin

# SingleFault class
class SingleFault(FaultsBin):
  """
  Python faults container with one material.

  Factory: faults_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FaultsBin.Inventory):
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

    from Fault import Fault
    fault = pyre.inventory.facility("fault", family="fault", factory=Fault)
    fault.meta['tip'] = "Fault in problem."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fault"):
    """
    Constructor.
    """
    FaultsBin.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set attributes from inventory.
    """
    FaultsBin._configure(self)
    self.faults = [self.inventory.fault]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def faults_bin():
  """
  Factory associated with SingleFault.
  """
  return SingleFault()


# End of file 
