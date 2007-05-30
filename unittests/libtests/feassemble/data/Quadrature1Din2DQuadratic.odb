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

## @file unittests/libtests/feassemble/data/Quadrature1Din2DQuadratic.odb
##
## @brief Python container holding quadrature information for a 1-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from Quadrature1DQuadratic import *

# ----------------------------------------------------------------------

# Quadrature1Din2DQuadratic class
class Quadrature1Din2DQuadratic(Quadrature1DQuadratic):
  """
  Python container holding quadrature information for a 1-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature1din2dquadratic"):
    """
    Constructor.
    """
    Quadrature1DQuadratic.__init__(self, name)
    
    self.spaceDim = 2
    return


# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature1DQuadratic.
  """
  return Quadrature1Din2DQuadratic()


# End of file 
