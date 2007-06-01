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

## @file pylith/materials/ElasticPlaneStress.py
##
## @brief Python object implementing 2-D isotropic linear elastic
## material for plane stress.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial

# ElasticPlaneStress class
class ElasticPlaneStress(ElasticMaterial):
  """
  Python object implementing 2-D isotropic linear elastic material for
  plane stress.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticplanestress"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    import pylith.materials.materials as bindings
    self.cppHandle = bindings.ElasticPlaneStress()
    self.dimension = self.cppHandle.dimension
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticPlaneStress.
  """
  return ElasticPlaneStress()


# End of file 
