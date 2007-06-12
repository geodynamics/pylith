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

## @file pylith/bc/BCSixSides.py
##
## @brief Python boundary conditions container for a 3-D box.
##
## Boundary conditions can be applied to any of the six faces.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# BCSixSides class
class BCSixSides(ObjectBin):
  """
  Python boundary conditions container for a 3-D box.

  Boundary conditions can be applied to any of the six faces.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
    """
    Python object for managing BCSixSides facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BCSixSides facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b x_pos Boundary condition on +x face of 3-D box.
    ## @li \b x_neg Boundary condition on -x face of 3-D box.
    ## @li \b y_pos Boundary condition on +y face of 3-D box.
    ## @li \b y_neg Boundary condition on -y face of 3-D box.
    ## @li \b z_pos Boundary condition on +z face of 3-D box.
    ## @li \b z_neg Boundary condition on -z face of 3-D box.

    import pyre.inventory
    
    from pylith.bc.Dirichlet import Dirichlet

    xPos = pyre.inventory.facility("x_pos", family="boundary_condition",
                                   factory=Dirichlet)
    xPos.meta['tip'] = "Boundary condition on +x face of 3-D box."

    xNeg = pyre.inventory.facility("x_neg", family="boundary_condition",
                                   factory=Dirichlet)
    xNeg.meta['tip'] = "Boundary condition on -x face of 3-D box."

    yPos = pyre.inventory.facility("y_pos", family="boundary_condition",
                                   factory=Dirichlet)
    yPos.meta['tip'] = "Boundary condition on +y face of 3-D box."

    yNeg = pyre.inventory.facility("y_neg", family="boundary_condition",
                                   factory=Dirichlet)
    yNeg.meta['tip'] = "Boundary condition on -y face of 3-D box."

    zPos = pyre.inventory.facility("z_pos", family="boundary_condition",
                                   factory=Dirichlet)
    zPos.meta['tip'] = "Boundary condition on +z face of 3-D box."

    zNeg = pyre.inventory.facility("z_neg", family="boundary_condition",
                                   factory=Dirichlet)
    zNeg.meta['tip'] = "Boundary condition on -z face of 3-D box."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bcsixsides"):
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
                self.inventory.yNeg,
                self.inventory.zPos,
                self.inventory.zNeg]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with BCSixSides.
  """
  return BCSixSides()


# End of file 
