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

## @file unittests/libtests/feassemble/data/Solution2DLinear.py

## @brief Python container holding solution information for 2-D mesh
## with linear cells used in testing finite-element C++ routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Solution2DLinear class
class Solution2DLinear(Component):
  """
  Python container holding solution information for 2-D mesh with
  linear cells used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="solution2dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution")
    
    # Input fields
    self.dt = 0.01
    self.fieldTpdt = numpy.array([ 1.9, -0.9,
                                   1.4,  1.5,
                                   0.5, -0.9 ], dtype=numpy.float64)
    self.fieldT = numpy.array([ 1.6, -0.8,
                                0.9,  0.7,
                               -0.2, -1.1 ], dtype=numpy.float64)
    self.fieldTmdt = numpy.array([ 0.8,  0.1,
                                   0.5,  0.3,
                                  -0.1, -0.6 ], dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def solution():
  """
  Factory for Solution2DLinear.
  """
  return Solution2DLinear()


# End of file 
