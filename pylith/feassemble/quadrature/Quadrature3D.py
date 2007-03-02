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

## @file pylith/feassemble/Qudrature3D.py
##
## @brief Python object implementing 3-D integration in 3-D space
## using numerical quadrature.
##
## Factory: quadrature

from Quadrature import Quadrature

# Quadrature3D class
class Quadrature3D(Quadrature):
  """
  Python object for integrating over 3-D finite-elements in a 3-D
  domain using quadrature.

  Factory: quadrature.
  """

  def __init__(self, name="quadrature3d"):
    """
    Constructor.
    """
    Quadrature.__init__(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature3D()
    self.spaceDim = 3
    self.cellDim = 3
    return


# FACTORIES ////////////////////////////////////////////////////////////

def quadrature():
  """
  Factory associated with Quadrature3D.
  """
  return Quadrature3D()


# End of file 
