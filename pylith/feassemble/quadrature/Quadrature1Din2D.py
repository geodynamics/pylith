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

## @file pylith/feassemble/quadrature/Qudrature1Din2D.py
##
## @brief Python object implementing 1-D integration in 2-D space
## using numerical quadrature.
##
## Factory: quadrature

from Quadrature import Quadrature

# Quadrature1Din2D class
class Quadrature1Din2D(Quadrature):
  """
  Python object for integrating over 1-D finite-elements in a 2-D
  domain using quadrature.

  Factory: quadrature.
  """

  def __init__(self, name="quadrature1din2d"):
    """
    Constructor.
    """
    Quadrature.__init__(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature1Din2D()
    self.spaceDim = 2
    self.cellDim = 1
    return


# FACTORIES ////////////////////////////////////////////////////////////

def quadrature():
  """
  Factory associated with Quadrature1Din2D.
  """
  return Quadrature1Din2D()


# End of file 
