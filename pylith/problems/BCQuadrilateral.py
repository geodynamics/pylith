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

## @file pylith/problems/BCQuadrilateral.py
##
## @brief Python boundary conditions container for a 2-D quadrilateral.
##
## Boundary conditions can be applied to any of the four edges.
##
## Factory: boundary_conditions

from BoundaryConditions import BoundaryConditions

# BCQuadrilateral class
class BCQuadrilateral(BoundaryConditions):
  """
  Python boundary conditions container for a 2-D quadrilateral.

  Boundary conditions can be applied to any of the four edges.

  Factory: boundary_conditions
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(BoundaryConditions.Inventory):
    """
    Python object for managing BCQuadrilateral facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BCQuadrilateral facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b x_pos Boundary condition on +x face of prism.
    ## @li \b x_neg Boundary condition on -x face of prism.
    ## @li \b y_pos Boundary condition on +y face of prism.
    ## @li \b y_neg Boundary condition on -y face of prism.

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
                                   factory=BoundaryCondition)
    yPos.meta['tip'] = "Boundary condition on +y face of prism."

    yNeg = pyre.inventory.facility("y_neg", family="boundary_condition",
                                   factory=Dirichlet)
    yNeg.meta['tip'] = "Boundary condition on -y face of prism."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bcquadrilateral"):
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
               self.inventory.yNeg]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def boundary_conditions():
  """
  Factory associated with BCQuadrilateral.
  """
  return BCQuadrilateral()


# End of file 
