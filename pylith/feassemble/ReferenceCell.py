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
    self.basis = None
    self.basisDeriv = None
    self.quadrature = None
    self.cellDim = None
    self.numCorners = None
    self.numQuadPts = None
    return

  def initialize(self):
    """
    Initialize reference finite-element cell.
    """
    return


# End of file 
