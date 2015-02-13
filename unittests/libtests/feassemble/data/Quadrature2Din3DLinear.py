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

## @file unittests/libtests/feassemble/data/Quadrature2Din3DLinear.odb
##
## @brief Python container holding quadrature information for a 2-D
## linear finite-element cell used in testing finite-element C++
## routines.

from Quadrature2DLinear import *

# ----------------------------------------------------------------------

# Quadrature2Din3DLinear class
class Quadrature2Din3DLinear(Quadrature2DLinear):
  """
  Python container holding quadrature information for a 2-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature2din3dlinear"):
    """
    Constructor.
    """
    Quadrature2DLinear.__init__(self, name)
    
    self.spaceDim = 3
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature2Din3DLinear.
  """
  return Quadrature2Din3DLinear()


# End of file 
