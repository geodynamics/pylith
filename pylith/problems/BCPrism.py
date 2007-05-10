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

## @file pylith/problems/BCPrism.py
##
## @brief Python boundary conditions container for a 3-D prism.
##
## Boundary conditions can be applied to any of the six faces.
##
## Factory: boundary_conditions

from BoundaryConditions import BoundaryConditions

# BCPrism class
class BCPrism(BoundaryConditions):
  """
  Python boundary conditions container for a 3-D prism.

  Boundary conditions can be applied to any of the six faces.

  Factory: boundary_conditions
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(BoundaryConditions.Inventory):
    """
    Python object for managing BCPrism facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BCPrism facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b x_pos Boundary condition on +x face of prism.
    ## @li \b x_neg Boundary condition on -x face of prism.
    ## @li \b y_pos Boundary condition on +y face of prism.
    ## @li \b y_neg Boundary condition on -y face of prism.
    ## @li \b z_pos Boundary condition on +z face of prism.
    ## @li \b z_neg Boundary condition on -z face of prism.

    import pyre.inventory
    
    from pylith.bc.Dirichlet import Dirichlet
    from pylith.bc.BoundaryCondition import BoundaryCondition

    xPos = pyre.inventory.facility("x_pos", family="boundary_condition",
                                   factory=Dirichlet)
    xPos.meta['tip'] = "Boundary condition on +x face of prism."

    xNeg = pyre.inventory.facility("x_neg", family="boundary_condition",
                                   factory=Dirichlet)
    xNeg.meta['tip'] = "Boundary condition on -x face of prism."

    yPos = pyre.inventory.facility("y_pos", family="boundary_condition",
                                   factory=Dirichlet)
    yPos.meta['tip'] = "Boundary condition on +y face of prism."

    yNeg = pyre.inventory.facility("y_neg", family="boundary_condition",
                                   factory=Dirichlet)
    yNeg.meta['tip'] = "Boundary condition on -y face of prism."

    zPos = pyre.inventory.facility("z_pos", family="boundary_condition",
                                   factory=BoundaryCondition)
    zPos.meta['tip'] = "Boundary condition on +z face of prism."

    zNeg = pyre.inventory.facility("z_neg", family="boundary_condition",
                                   factory=Dirichlet)
    zNeg.meta['tip'] = "Boundary condition on -z face of prism."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bcprism"):
    """
    Constructor.
    """
    BoundaryConditions.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set attributes from inventory.
    """
    BoundaryConditions._configure(self)
    self.bc = [self.inventory.xPos,
               self.inventory.xNeg,
               self.inventory.yPos,
               self.inventory.yNeg,
               self.inventory.zPos,
               self.inventory.zNeg]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def boundary_conditions():
  """
  Factory associated with BCPrism.
  """
  return BCPrism()


# End of file 
