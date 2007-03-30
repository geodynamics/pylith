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

from Material import Material

import pylith.materials.materials as bindings


# ElasticPlaneStress class
class ElasticPlaneStress(Material):
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
    Material.__init__(self, name)
    self.cppHandle = bindings.ElasticPlaneStress()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticPlaneStress.
  """
  return ElasticPlaneStress()


# End of file 
