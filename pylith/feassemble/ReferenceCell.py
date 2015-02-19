#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/feassemble/ReferenceCell.py
##
## @brief Python abstract base class for managing basis functions and
## quadrature rules of a reference finite-element cell.
##
## Factory: reference_cell.

from pylith.utils.PetscComponent import PetscComponent

# ReferenceCell class
class ReferenceCell(PetscComponent):
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
    PetscComponent.__init__(self, name, facility="reference_cell")

    self.geometry = None # Geometry of reference cell

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


  def initialize(self, spaceDim):
    """
    Initialize reference finite-element cell.
    """
    return


# End of file 
