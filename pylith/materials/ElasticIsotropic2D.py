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
## @brief Python object implementing 2-D isotropic linear elastic material.
##
## Factory: material.

from Material import Material

# import pylith.materials.materials as bindings


# ElasticIsotropic2D class
class ElasticIsotropic2D(Material):
  """
  Python object implementing 2-D isotropic linear elastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticisotropic2d"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    # :TODO: Need to create module for materials
    # self.cppHandle = bindings.ElasticIsotropic2D()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticIsotropic2D.
  """
  return ElasticIsotropic2D()


# End of file 
