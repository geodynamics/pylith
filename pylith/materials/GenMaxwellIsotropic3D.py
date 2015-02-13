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

## @file pylith/materials/GenMaxwellIsotropic3D.py
##
## @brief Python object implementing 3-D isotropic linear GenMaxwell
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import GenMaxwellIsotropic3D as ModuleGenMaxwellIsotropic3D

# GenMaxwellIsotropic3D class
class GenMaxwellIsotropic3D(ElasticMaterial, ModuleGenMaxwellIsotropic3D):
  """
  Python object implementing 3-D isotropic linear GenMaxwell viscoelastic
  material.

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
           {'info': ["mu", "lambda", "density", "stable_dt_implicit", "stable_dt_explicit",
                     "shear_ratio", "maxwell_time"],
            'data': ["total_strain", "stress",
                     "viscous_strain_1", 
                     "viscous_strain_2", 
                     "viscous_strain_3",
                     ]}}
    self._loggingPrefix = "MaGM3D "
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleGenMaxwellIsotropic3D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with GenMaxwellIsotropic3D.
  """
  return GenMaxwellIsotropic3D()


# End of file 
