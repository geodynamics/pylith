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

## @file pylith/materials/ElasticStrain1D.py
##
## @brief Python object implementing 1-D linear elastic material with
## axial strain.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial

# ElasticStrain1D class
class ElasticStrain1D(ElasticMaterial):
  """
  Python object implementing 1-D linear elastic material with axial strain.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticstrain1d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    import pylith.materials.materials as bindings
    self.cppHandle = bindings.ElasticStrain1D()
    self.dimension = self.cppHandle.dimension
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticStrain1D.
  """
  return ElasticStrain1D()


# End of file 
