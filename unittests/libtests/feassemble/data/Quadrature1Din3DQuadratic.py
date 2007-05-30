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

## @file unittests/libtests/feassemble/data/Quadrature1Din3DQuadratic.odb
##
## @brief Python container holding quadrature information for a 1-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from Quadrature1DQuadratic import *

# ----------------------------------------------------------------------

# Quadrature1Din3DQuadratic class
class Quadrature1Din3DQuadratic(Quadrature1DQuadratic):
  """
  Python container holding quadrature information for a 1-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature1din3dquadratic"):
    """
    Constructor.
    """
    Quadrature1DQuadratic.__init__(self, name)
    
    self.spaceDim = 3
    return


# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature1Din3DQuadratic.
  """
  return Quadrature1Din3DQuadratic()


# End of file 
