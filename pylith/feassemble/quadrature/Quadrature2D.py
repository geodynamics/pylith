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

## @file pylith/feassemble/quadrature/Qudrature2D.py
##
## @brief Python object implementing 2-D integration in 2-D space
## using numerical quadrature.
##
## Factory: quadrature

from Quadrature import Quadrature

# Quadrature2D class
class Quadrature2D(Quadrature):
  """
  Python object for integrating over 2-D finite-elements in a 2-D
  domain using quadrature.

  Factory: quadrature.
  """

  def __init__(self, name="quadrature2d"):
    """
    Constructor.
    """
    Quadrature.__init__(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature2D()
    self.spaceDim = 2
    self.cellDim = 2
    return


# FACTORIES ////////////////////////////////////////////////////////////

def quadrature():
  """
  Factory associated with Quadrature2D.
  """
  return Quadrature2D()


# End of file 
