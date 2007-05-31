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

## @file unittests/libtests/feassemble/data/Solution3DQuadratic.py

## @brief Python container holding solution information for 3-D mesh
## with quadratic cells used in testing finite-element C++ routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Solution3DQuadratic class
class Solution3DQuadratic(Component):
  """
  Python container holding solution information for 3-D mesh with
  quadratic cells used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="solution3dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution")
    
    # Input fields
    self.dt = 0.01
    self.fieldTpdt = numpy.array([ +0.3, -0.4, -0.4,
                                   -0.6, +0.8, +0.2,
                                   +0.5, +0.5, +0.7,
                                   -0.7, -0.5, -0.7,
                                   -0.6, -0.3, +0.8,
                                   -0.4, -0.8, -0.5,
                                   +0.7, +0.8, -0.5,
                                   -0.5, -0.5, -0.7,
                                   -0.3, -0.9, +0.8,
                                   -0.1, +0.5, -0.9], dtype=numpy.float64)
    self.fieldT = numpy.array([ +0.1, -0.2, -0.6,
                                -0.3, +0.4, +0.9,
                                +0.6, +0.8, +0.5,
                                -0.8, -0.6, -0.8,
                                -0.0, -0.2, +0.6,
                                -0.4, -0.7, -0.2,
                                +0.7, +0.6, -0.1,
                                -0.4, -0.3, -0.3,
                                -0.7, -0.6, +0.1,
                                -0.9, +0.3, -0.8], dtype=numpy.float64)
    self.fieldTmdt = numpy.array([ +0.2, -0.3, -0.4,
                                   -0.6, +0.2, +0.3,
                                   +0.5, +0.2, +0.5,
                                   -0.3, -0.6, -0.3,
                                   -0.5, -0.9, +0.4,
                                   -0.3, -0.6, -0.8,
                                   +0.9, +0.5, -0.2,
                                   -0.7, -0.3, -0.3,
                                   -0.5, -0.8, +0.4,
                                   -0.4, +0.5, -0.8], dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def solution():
  """
  Factory for Solution3DQuadratic.
  """
  return Solution3DQuadratic()


# End of file 
