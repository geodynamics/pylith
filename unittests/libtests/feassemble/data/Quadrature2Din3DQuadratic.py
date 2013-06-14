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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Quadrature2Din3DQuadratic.odb
##
## @brief Python container holding quadrature information for a 2-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from Quadrature2DQuadratic import *

# ----------------------------------------------------------------------

# Quadrature2Din3DQuadratic class
class Quadrature2Din3DQuadratic(Quadrature2DQuadratic):
  """
  Python container holding quadrature information for a 1-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature2din3dquadratic"):
    """
    Constructor.
    """
    Quadrature2DQuadratic.__init__(self, name)
    
    self.spaceDim = 3
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature2Din3DQuadratic.
  """
  return Quadrature2Din3DQuadratic()


# End of file 
