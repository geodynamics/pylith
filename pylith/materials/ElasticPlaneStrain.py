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

from ElasticMaterial import ElasticMaterial

# ElasticPlaneStrain class
class ElasticPlaneStrain(ElasticMaterial):
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
    ElasticMaterial.__init__(self, name)
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.materials.materials as bindings
      self.cppHandle = bindings.ElasticPlaneStrain()
      self.dimension = self.cppHandle.dimension
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticPlaneStrain.
  """
  return ElasticPlaneStrain()


# End of file 
