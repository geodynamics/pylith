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

## @file pylith/feassemble/FIATCell.py

## @brief Python object for managing basis functions and quadrature
## rules of a reference finite-element cell using FIAT.

from ReferenceCell import ReferenceCell

import numpy

# FIATCell class
class FIATCell(ReferenceCell):
  """
  Python object for managing basis functions and quadrature rules of a
  reference finite-element cell using FIAT.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fiatcell"):
    """
    Constructor.
    """
    ReferenceCell.__init__(self, name)
    return

  def initialize(self):
    """
    Initialize reference finite-element cell.
    """
    quadrature = self._setupQuadrature()
    basisFns = self._setupBasisFns()
    
    # Evaluate basis functions at quadrature points
    quadpts = quadrature.get_points()
    self.basis = numpy.array(basisFns.tabulate(quadpts)).transpose()

    # Evaluate derivatives of basis functions at quadrature points
    # ADD STUFF HERE

    self._info.line("Basis:")
    self._info.line(self.basis)
    #print "Basis derivatives:"
    #print self.basisDeriv
    self._info.line("Quad pts:")
    self._info.line(quadrature.get_points())
    self._info.line("Quad wts:")
    self._info.line(quadrature.get_weights())
    self._info.log()

    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """
    raise NotImplementedError()
    return


  def _setupBasisFns(self):
    """
    Setup basis functions for reference cell.
    """
    raise NotImplementedError()
    return


# End of file 
