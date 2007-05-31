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
    self.lameMu = 3.0e+10
    self.lameLambda = self.lameMu
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def material():
  """
  Factory for ElasticPlaneStrain.
  """
  return ElasticPlaneStrain()


# End of file 
