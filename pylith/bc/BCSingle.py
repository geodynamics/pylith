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

## @file pylith/bc/BCSingle.py
##
## @brief Python boundary conditions container for one boundary condition.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# BCSingle class
class BCSingle(ObjectBin):
  """
  Python boundary conditions container for one boundary condition.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
    """
    Python object for managing BCSingle facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BCSingle facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b bc Boundary condition.

    import pyre.inventory
    
    from pylith.bc.Dirichlet import Dirichlet

    bc = pyre.inventory.facility("bc", family="boundary_condition",
                                 factory=Dirichlet)
    bc.meta['tip'] = "Boundary condition."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bcsingle"):
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
    self.bin = [self.inventory.bc]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with BCSingle.
  """
  return BCSingle()


# End of file 
