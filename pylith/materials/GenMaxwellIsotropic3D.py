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

## @file pylith/materials/GenMaxwellIsotropic3D.py
##
## @brief Python object implementing 3-D generalized isotropic linear Maxwell viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial

# GenMaxwellIsotropic3D class
class GenMaxwellIsotropic3D(ElasticMaterial):
  """
  Python object implementing 3-D generalized isotropic linear Maxwell viscoelastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="genmaxwellisotropic3d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["mu", "lambda", "density", "shear-ratio", "Maxwell-time"],
            'data': ["total-strain", "viscous-strain", "stress"]}}
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.materials.materials as bindings
      self.cppHandle = bindings.GenMaxwellIsotropic3D()
      self.dimension = self.cppHandle.dimension
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with GenMaxwellIsotropic3D.
  """
  return GenMaxwellIsotropic3D()


# End of file 
