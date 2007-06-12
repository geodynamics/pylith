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

## @file pylith/bc/BCFourSides.py
##
## @brief Python boundary conditions container for a 2-D quadrilateral.
##
## Boundary conditions can be applied to any of the four edges.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# BCFourSides class
class BCFourSides(ObjectBin):
  """
  Python boundary conditions container for a 2-D quadrilateral.

  Boundary conditions can be applied to any of the four edges.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
    """
    Python object for managing BCFourSides facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BCFourSides facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b x_pos Boundary condition on +x face of 2-D box.
    ## @li \b x_neg Boundary condition on -x face of 2-D box.
    ## @li \b y_pos Boundary condition on +y face of 2-D box.
    ## @li \b y_neg Boundary condition on -y face of 2-D box.

    import pyre.inventory
    
    from pylith.bc.Dirichlet import Dirichlet

    xPos = pyre.inventory.facility("x_pos", family="boundary_condition",
                                   factory=Dirichlet)
    xPos.meta['tip'] = "Boundary condition on +x face of 2-D box."

    xNeg = pyre.inventory.facility("x_neg", family="boundary_condition",
                                   factory=Dirichlet)
    xNeg.meta['tip'] = "Boundary condition on -x face of 2-D box."

    yPos = pyre.inventory.facility("y_pos", family="boundary_condition",
                                   factory=Dirichlet)
    yPos.meta['tip'] = "Boundary condition on +y face of 2-D box."

    yNeg = pyre.inventory.facility("y_neg", family="boundary_condition",
                                   factory=Dirichlet)
    yNeg.meta['tip'] = "Boundary condition on -y face of 2-D box."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bcfoursides"):
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
    self.bin = [self.inventory.xPos,
                self.inventory.xNeg,
                self.inventory.yPos,
                self.inventory.yNeg]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with BCFourSides.
  """
  return BCFourSides()


# End of file 
