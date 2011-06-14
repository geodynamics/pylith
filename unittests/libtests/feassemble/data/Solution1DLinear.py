#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Solution1DLinear.py

## @brief Python container holding solution information for 1-D mesh
## with linear cells used in testing finite-element C++ routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------

# Solution1DLinear class
class Solution1DLinear(Component):
  """
  Python container holding solution information for 1-D mesh with
  linear cells used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="solution1dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution")
    
    # Input fields
    self.dt = 0.01
    self.fieldTIncr = numpy.array( [[1.2], [1.7]], dtype=numpy.float64)
    self.fieldT = numpy.array( [[1.1], [1.5]], dtype=numpy.float64)
    self.fieldTmdt = numpy.array( [[1.0], [1.3]], dtype=numpy.float64)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def solution():
  """
  Factory for Solution1DLinear.
  """
  return Solution1DLinear()


# End of file 
