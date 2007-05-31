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

## @file unittests/libtests/feassemble/data/ElasticIsotropic3D.py

## @brief Python container holding material information for elastic
## isotropic 3-D material used in testing finite-element C++ routines.

from pyre.components.Component import Component

# ----------------------------------------------------------------------

# ElasticIsotropic3D class
class ElasticIsotropic3D(Component):
  """
  Python container holding material information for elastic isotropic 3-D
  material used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticisotropic3d"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="material")
    
    self.dimension = 3
    self.type = "ElasticIsotropic3D"
    self.dbFilename = "data/elasticisotropic3d.spatialdb"
    self.id = 0
    self.label = "elastic isotropic 3-D"
    self.density = 2500.0
    self.lameMu = 3.0e+10
    self.lameLambda = self.lameMu
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def material():
  """
  Factory for ElasticIsotropic3D.
  """
  return ElasticIsotropic3D()


# End of file 
