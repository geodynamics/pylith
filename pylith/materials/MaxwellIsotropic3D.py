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

## @file pylith/materials/MaxwellIsotropic3D.py
##
## @brief Python object implementing 3-D isotropic linear Maxwell viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial

# MaxwellIsotropic3D class
class MaxwellIsotropic3D(ElasticMaterial):
  """
  Python object implementing 3-D isotropic linear Maxwell viscoelastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="maxwellisotropic3d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    import pylith.materials.materials as bindings
    self.cppHandle = bindings.MaxwellIsotropic3D()
    self.dimension = self.cppHandle.dimension
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with MaxwellIsotropic3D.
  """
  return MaxwellIsotropic3D()


# End of file 
