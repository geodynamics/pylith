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
## @brief Python object implementing 1-D isotropic linear elastic material.
##
## Factory: material.

from Material import Material

# import pylith.materials.materials as bindings


# ElasticIsotropic1D class
class ElasticIsotropic3D(Material):
  """
  Python object implementing 1-D isotropic linear elastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticisotropic1d"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    # :TODO: Need to create module for materials
    # self.cppHandle = bindings.ElasticIsotropic1D()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticIsotropic1D.
  """
  return ElasticIsotropic1D()


# End of file 
