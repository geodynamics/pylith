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

## @file pylith/materials/ElasticPlaneStrain.py
##
## @brief Python object implementing 1-D isotropic linear elastic
## material for plane strain.
##
## Factory: material.

from Material import Material

import pylith.materials.materials as bindings


# ElasticPlaneStrain class
class ElasticPlaneStrain(Material):
  """
  Python object implementing 2-D isotropic linear elastic material for
  plane strain.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticplanestrain"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    self.cppHandle = bindings.ElasticPlaneStrain()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticPlaneStrain.
  """
  return ElasticPlaneStrain()


# End of file 
