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

## @file unittests/libtests/feassemble/data/Solution3DLinear.py

## @brief Python container holding solution information for 3-D mesh
## with linear cells used in testing finite-element C++ routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Solution3DLinear class
class Solution3DLinear(Component):
  """
  Python container holding solution information for 3-D mesh with
  linear cells used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="solution3dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution")
    
    # Input fields
    self.dt = 0.01
    self.fieldTpdt = numpy.array([ +0.3, +0.2, -0.5,
                                   -0.3, -0.4, -0.6,
                                   +0.2, +0.6, +0.3,
                                   -0.6, -0.1, -0.3], dtype=numpy.float64)
    self.fieldT = numpy.array([ +0.8, +0.1, -0.6,
                                -0.1, -0.2, -0.5,
                                +0.1, +0.7, +0.2,
                                -0.5, -0.0, -0.2], dtype=numpy.float64)
    self.fieldTmdt = numpy.array([ +0.1, +0.1, -0.3,
                                   -0.2, -0.1, -0.5,
                                   +0.2, +0.4, +0.1,
                                   -0.4, -0.1, -0.1], dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def solution():
  """
  Factory for Solution3DLinear.
  """
  return Solution3DLinear()


# End of file 
