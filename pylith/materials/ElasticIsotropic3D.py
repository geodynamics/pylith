#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/ElasticIsotropic1D.py
##
## @brief Python object implementing 3-D isotropic linear elastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import ElasticIsotropic3D as ModuleElasticIsotropic3D

# ElasticIsotropic3D class
class ElasticIsotropic3D(ElasticMaterial, ModuleElasticIsotropic3D):
  """
  Python object implementing 3-D isotropic linear elastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticisotropic3d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["mu", "lambda", "density", "stable_dt_implicit", "stable_dt_explicit",],
            'data': ["total_strain", "stress"]}}
    self._loggingPrefix = "MaEl3D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleElasticIsotropic3D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticIsotropic3D.
  """
  return ElasticIsotropic3D()


# End of file 
