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

from Material import Material

import pylith.materials.materials as bindings


# ElasticStrain1D class
class ElasticStrain1D(Material):
  """
  Python object implementing 1-D linear elastic material with axial strain.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticstrain1d"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    self.cppHandle = bindings.ElasticStrain1D()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticStrain1D.
  """
  return ElasticStrain1D()


# End of file 
