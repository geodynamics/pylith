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

from ElasticMaterial import ElasticMaterial

# ElasticStress1D class
class ElasticStress1D(ElasticMaterial):
  """
  Python object implementing 1-D linear elastic material with axial stress.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticstress1d"):
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
      self.cppHandle = bindings.ElasticStress1D()
      self.dimension = self.cppHandle.dimension
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticStress1D.
  """
  return ElasticStress1D()


# End of file 
