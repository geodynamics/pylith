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

## @file unittests/libtests/feassemble/data/ElasticPlaneStrain.py
##
## @brief Python container holding material information for elastic
## strain 2-D material used in testing finite-element C++ routines.

from pyre.components.Component import Component

# ----------------------------------------------------------------------

# ElasticPlaneStrain class
class ElasticPlaneStrain(Component):
  """
  Python container holding material information for elastic strain 2-D
  material used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticplanestrain"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="material")
    
    self.dimension = 2
    self.type = "ElasticPlaneStrain"
    self.dbFilename = "data/elasticplanestrain.spatialdb"
    self.id = 0
    self.label = "elastic strain 2-D"
    self.density = 2500.0
    vs = 3500.0
    vp = 6000.0
    self.lameMu = self.density*vs*vs
    self.lameLambda = self.density*vp*vp - 2.0*self.lameMu
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def material():
  """
  Factory for ElasticPlaneStrain.
  """
  return ElasticPlaneStrain()


# End of file 
