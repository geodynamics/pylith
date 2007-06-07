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

## @file pylith/feassemble/quadrature/Qudrature2Din3D.py
##
## @brief Python object implementing 2-D integration in 3-D space
## using numerical quadrature.
##
## Factory: quadrature

from Quadrature import Quadrature

# Quadrature2Din3D class
class Quadrature2Din3D(Quadrature):
  """
  Python object for integrating over 2-D finite-elements in a 3-D
  domain using quadrature.

  Factory: quadrature.
  """

  def __init__(self, name="quadrature2din3d"):
    """
    Constructor.
    """
    Quadrature.__init__(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature2Din3D()
    self.spaceDim = 3
    self.cellDim = 2
    return


# FACTORIES ////////////////////////////////////////////////////////////

def quadrature():
  """
  Factory associated with Quadrature2Din3D.
  """
  return Quadrature2Din3D()


# End of file 
