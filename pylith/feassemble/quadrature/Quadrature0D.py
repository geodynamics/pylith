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

## @file pylith/feassemble/quadrature/Qudrature0D.py
##
## @brief Python object implementing 0-D integration in 1-D space
## using numerical quadrature.
##
## Factory: quadrature

from Quadrature import Quadrature

# Quadrature0D class
class Quadrature0D(Quadrature):
  """
  Python object for integrating over 0-D finite-elements in a 1-D
  domain using quadrature.

  Factory: quadrature.
  """

  def __init__(self, name="quadrature0d"):
    """
    Constructor.
    """
    Quadrature.__init__(self, name)
    self.spaceDim = 1
    self.cellDim = 0
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.feassemble.feassemble as bindings
      self.cppHandle = bindings.Quadrature0D()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def quadrature():
  """
  Factory associated with Quadrature0D.
  """
  return Quadrature0D()


# End of file 
