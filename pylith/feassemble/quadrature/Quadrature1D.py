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

## @file pylith/feassemble/quadrature/Qudrature1D.py
##
## @brief Python object implementing 1-D integration in 1-D space
## using numerical quadrature.
##
## Factory: quadrature

from Quadrature import Quadrature

# Quadrature1D class
class Quadrature1D(Quadrature):
  """
  Python object for integrating over 1-D finite-elements in a 1-D
  domain using quadrature.

  Factory: quadrature.
  """

  def __init__(self, name="quadrature1d"):
    """
    Constructor.
    """
    Quadrature.__init__(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature1D()
    self.spaceDim = 1
    self.cellDim = 1
    return


# FACTORIES ////////////////////////////////////////////////////////////

def quadrature():
  """
  Factory associated with Quadrature1D.
  """
  return Quadrature1D()


# End of file 
