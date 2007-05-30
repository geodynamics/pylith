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

## @file unittests/libtests/feassemble/data/MaterialElasticStrain1D.py

## @brief Python container holding material information for elastic
## strain 1-D material used in testing finite-element C++ routines.

from pyre.components.Component import Component

# ----------------------------------------------------------------------

# MaterialElasticStrain1D class
class MaterialElasticStrain1D(Component):
  """
  Python container holding material information for elastic strain 1-D
  material used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="materialelasticstrain1d"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="material")
    
    self.dimension = 1
    self.type = "ElasticStrain1D"
    self.dbFilename = "data/elasticstrain1d.spatialdb"
    self.id = 0
    self.label = "elastic strain 1-D"
    self.density = 2500.0
    self.lameMu = 3.0e+10
    self.lameLambda = self.lameMu
    return
  

# FACTORIES ////////////////////////////////////////////////////////////
def material():
  """
  Factory for MaterialElasticStrain1D.
  """
  return MaterialElasticStrain1D()


# End of file 
