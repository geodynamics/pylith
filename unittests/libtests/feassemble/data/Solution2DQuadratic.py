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

## @file unittests/libtests/feassemble/data/Solution2DQuadratic.py

## @brief Python container holding solution information for 2-D mesh
## with quadratic cells used in testing finite-element C++ routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Solution2DQuadratic class
class Solution2DQuadratic(Component):
  """
  Python container holding solution information for 2-D mesh with
  quadratic cells used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="solution2dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution")
    
    # Input fields
    self.dt = 0.01
    self.fieldTpdt = numpy.array([ -0.4, -0.6,
                                   +0.7, +0.8,
                                   +0.0, +0.2,
                                   -0.5, -0.4,
                                   +0.3, +0.9,
                                   -0.3, -0.9 ], dtype=numpy.float64)
    self.fieldT = numpy.array([ -0.3, -0.4,
                                +0.5, +0.6,
                                +0.0, +0.1,
                                -0.2, -0.3,
                                +0.2, +0.3,
                                -0.1, -0.2 ], dtype=numpy.float64)
    self.fieldTmdt = numpy.array([ -0.2, -0.3,
                                   +0.3, +0.4,
                                   +0.0, +0.1,
                                   -0.3, -0.2,
                                   +0.1, +0.4,
                                   -0.2, -0.6 ], dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def solution():
  """
  Factory for Solution2DQuadratic.
  """
  return Solution2DQuadratic()


# End of file 
