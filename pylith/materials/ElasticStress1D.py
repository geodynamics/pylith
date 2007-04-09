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

## @file pylith/materials/ElasticStress1D.py
##
## @brief Python object implementing 1-D linear elastic material with
## axial stress.
##
## Factory: material.

from Material import Material

import pylith.materials.materials as bindings


# ElasticStress1D class
class ElasticStress1D(Material):
  """
  Python object implementing 1-D linear elastic material with axial stress.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticstress1d"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    self.cppHandle = bindings.ElasticStress1D()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticStress1D.
  """
  return ElasticStress1D()


# End of file 
