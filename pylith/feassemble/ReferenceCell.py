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

## @brief Python abstract base class for managing basis functions and
## quadrature rules of a reference finite-element cell.

from pyre.components.Component import Component

# ReferenceCell class
class ReferenceCell(Component):
  """
  Python abstract base class for managing basis functions and
  quadrature rules of a reference finite-element cell.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="referencecell"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="referencecell")

    self.basis = None # numpy array w/basis fns at quad pts
    self.basisDeriv = None # numpy array w/basis fn derivs at quad pts
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
