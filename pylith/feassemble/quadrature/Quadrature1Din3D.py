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

## @file pylith/feassemble/quadrature/Qudrature1Din3D.py
##
## @brief Python object implementing 1-D integration in 3-D space
## using numerical quadrature.
##
## Factory: quadrature

from Quadrature import Quadrature

# Quadrature1Din3D class
class Quadrature1Din3D(Quadrature):
  """
  Python object for integrating over 1-D finite-elements in a 3-D
  domain using quadrature.

  Factory: quadrature.
  """

  def __init__(self, name="quadrature1din3d"):
    """
    Constructor.
    """
    Quadrature.__init__(self, name)
    self.spaceDim = 3
    self.cellDim = 1
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.feassemble.feassemble as bindings
      self.cppHandle = bindings.Quadrature1Din3D()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def quadrature():
  """
  Factory associated with Quadrature1D.
  """
  return Quadrature1D()


# End of file 
