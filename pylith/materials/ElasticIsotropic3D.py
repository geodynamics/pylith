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

## @file pylith/materials/ElasticIsotropic1D.py
##
## @brief Python object implementing 3-D isotropic linear elastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial

# ElasticIsotropic3D class
class ElasticIsotropic3D(ElasticMaterial):
  """
  Python object implementing 3-D isotropic linear elastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticisotropic3d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    import pylith.materials.materials as bindings
    self.cppHandle = bindings.ElasticIsotropic3D()
    self.dimension = self.cppHandle.dimension
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticIsotropic3D.
  """
  return ElasticIsotropic3D()


# End of file 
