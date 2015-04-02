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

## @file pylith/materials/MaxwellIsotropic3D.py
##
## @brief Python object implementing 3-D isotropic linear Maxwell
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import MaxwellIsotropic3D as ModuleMaxwellIsotropic3D

# MaxwellIsotropic3D class
class MaxwellIsotropic3D(ElasticMaterial, ModuleMaxwellIsotropic3D):
  """
  Python object implementing 3-D isotropic linear Maxwell viscoelastic
  material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="maxwellisotropic3d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["mu", "lambda", "density", "stable_dt_implicit", "stable_dt_explicit", "maxwell_time"],
            'data': ["total_strain", "viscous_strain", "stress"]}}
    self._loggingPrefix = "MaMx3D "
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleMaxwellIsotropic3D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with MaxwellIsotropic3D.
  """
  return MaxwellIsotropic3D()


# End of file 
