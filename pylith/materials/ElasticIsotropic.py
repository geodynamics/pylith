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

## @file pylith/materials/ElasticIsotropic.py

## @brief Python objects implementing isotropic linear elastic materials.

from Material import Material

# import pylith.materials.materials as bindings


# ======================================================================
# ElasticIsotropic1D class

## @brief Python object implementing 1-D isotropic linear elastic material.
class ElasticIsotropic3D(Material):
  """
  Python object implementing 1-D isotropic linear elastic material.
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


# ======================================================================
# ElasticIsotropic2D class

## @brief Python object implementing 2-D isotropic linear elastic material.
class ElasticIsotropic2D(Material):
  """
  Python object implementing 2-D isotropic linear elastic material.
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


# ======================================================================
# ElasticIsotropic3D class

## @brief Python object implementing 3-D isotropic linear elastic material.
class ElasticIsotropic3D(Material):
  """
  Python object implementing 3-D isotropic linear elastic material.
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


# End of file 
