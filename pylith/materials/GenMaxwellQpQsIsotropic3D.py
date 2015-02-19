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

## @file pylith/materials/GenMaxwellQpQsIsotropic3D.py
##
## @brief Python object implementing 3-D isotropic linear GenMaxwellQpQs
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import GenMaxwellQpQsIsotropic3D as ModuleGenMaxwellQpQsIsotropic3D

# GenMaxwellQpQsIsotropic3D class
class GenMaxwellQpQsIsotropic3D(ElasticMaterial, ModuleGenMaxwellQpQsIsotropic3D):
  """
  Python object implementing 3-D isotropic linear GenMaxwellQpQs viscoelastic
  material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="genmaxwellqpqsisotropic3d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["mu", "k", "density", "stable_dt_implicit", "stable_dt_explicit",
                     "shear_ratio", 
                     "bulk_ratio",
                     "maxwell_time_shear",
                     "maxwell_time_bulk"],
            'data': ["total_strain", "stress",
                     "viscous_deviatoric_strain", 
                     "viscous_mean_strain", 
                     ]}}
    self._loggingPrefix = "MaGQ3D "
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleGenMaxwellQpQsIsotropic3D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with GenMaxwellQpQsIsotropic3D.
  """
  return GenMaxwellQpQsIsotropic3D()


# End of file 
