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

## @file unittests/libtests/feassemble/data/Solution1DQuadratic.py

## @brief Python container holding solution information for 1-D mesh
## with quadratic cells used in testing finite-element C++ routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Solution1DQuadratic class
class Solution1DQuadratic(Component):
  """
  Python container holding solution information for 1-D mesh with
  quadratic cells used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="solution1dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution")
    
    # Input fields
    self.dt = 0.01
    self.fieldTpdt = numpy.array([ [1.2], [0.0], [1.7] ], dtype=numpy.float64)
    self.fieldT = numpy.array([ [1.1], [0.1], [1.5] ], dtype=numpy.float64)
    self.fieldTmdt = numpy.array([ [1.0], [0.1], [1.3] ], dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def solution():
  """
  Factory for Solution1DQuadratic.
  """
  return Solution1DQuadratic()


# End of file 
