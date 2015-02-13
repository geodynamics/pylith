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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Quadrature1Din2DLinear.odb
##
## @brief Python container holding quadrature information for a 1-D
## linear finite-element cell used in testing finite-element C++
## routines.

from Quadrature1DLinear import *

# ----------------------------------------------------------------------

# Quadrature1Din2DLinear class
class Quadrature1Din2DLinear(Quadrature1DLinear):
  """
  Python container holding quadrature information for a 1-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature1din2dlinear"):
    """
    Constructor.
    """
    Quadrature1DLinear.__init__(self, name)
    
    self.spaceDim = 2
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature1Din2DLinear.
  """
  return Quadrature1Din2DLinear()


# End of file 
