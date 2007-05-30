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

## @file unittests/libtests/feassemble/data/Quadrature1Din3DLinear.odb
##
## @brief Python container holding quadrature information for a 1-D
## linear finite-element cell used in testing finite-element C++
## routines.

from Quadrature1DLinear import *

import numpy

# ----------------------------------------------------------------------

# Quadrature1Din3DLinear class
class Quadrature1Din3DLinear(Quadrature1DLinear):
  """
  Python container holding quadrature information for a 1-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature1din3dlinear"):
    """
    Constructor.
    """
    Quadrature1DLinear.__init__(self, name)
    
    self.spaceDim = 3
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature1Din3DLinear.
  """
  return Quadrature1Din3DLinear()


# End of file 
