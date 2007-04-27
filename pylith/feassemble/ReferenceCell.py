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

## @file pylith/feassemble/ReferenceCell.py
##
## @brief Python abstract base class for managing basis functions and
## quadrature rules of a reference finite-element cell.
##
## Factory: reference_cell.

from pyre.components.Component import Component

# ReferenceCell class
class ReferenceCell(Component):
  """
  Python abstract base class for managing basis functions and
  quadrature rules of a reference finite-element cell.

  Factory: reference_cell.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="referencecell"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="reference_cell")

    self.basisVert = None # numpy array w/basis fns at vertices
    self.basisDerivVert = None # numpy array w/basis fn derivs at vertices
    self.basisQuad = None # numpy array w/basis fns at quad pts
    self.basisDerivQuad = None # numpy array w/basis fn derivs at quad pts

    self.quadPts = None # numpy array w/coordinates of quad pts
    self.quadWts = None # numpy array w/wts of quad pts

    self.cellDim = None # dimension of reference cell
    self.numCorners = None # number of vertices in reference cell
    self.numQuadPts = None # number of quadrature points
    return


  def initialize(self):
    """
    Initialize reference finite-element cell.
    """
    return


# End of file 
