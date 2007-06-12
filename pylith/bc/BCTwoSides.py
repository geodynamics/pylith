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

## @file pylith/bc/BCTwoSides.py
##
## @brief Python boundary conditions container for a 1-D box.
##
## Boundary conditions can be applied to any two points.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# BCTwoSides class
class BCTwoSides(ObjectBin):
  """
  Python boundary conditions container for a 1-D box.

  Boundary conditions can be applied to any two points.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
    """
    Python object for managing BCTwoSides facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BCTwoSides facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b pos Boundary condition on +x face of 1-D box.
    ## @li \b neg Boundary condition on -x face of 1-D box.

    import pyre.inventory
    
    from pylith.bc.Dirichlet import Dirichlet

    pos = pyre.inventory.facility("pos", family="boundary_condition",
                                  factory=Dirichlet)
    pos.meta['tip'] = "Boundary condition on positive face of 1-D box."
    
    neg = pyre.inventory.facility("neg", family="boundary_condition",
                                  factory=Dirichlet)
    neg.meta['tip'] = "Boundary condition on negative face of 1-D box."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bctwosides"):
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
    self.bin = [self.inventory.pos,
                self.inventory.neg]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with BCTwoSides.
  """
  return BCTwoSides()


# End of file 
