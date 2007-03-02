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

from Material import Material

# import pylith.materials.materials as bindings


# ElasticIsotropic3D class
class ElasticIsotropic3D(Material):
  """
  Python object implementing 3-D isotropic linear elastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticisotropic3d"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    # :TODO: Need to create module for materials
    # self.cppHandle = bindings.ElasticIsotropic3D()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticIsotropic3D.
  """
  return ElasticIsotropic3D()


# End of file 
