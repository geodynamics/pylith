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

## @file pylith/materials/GenMaxwellPlaneStrain.py
##
## @brief Python object implementing 3-D isotropic linear GenMaxwell
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import GenMaxwellPlaneStrain as ModuleGenMaxwellPlaneStrain

# GenMaxwellPlaneStrain class
class GenMaxwellPlaneStrain(ElasticMaterial, ModuleGenMaxwellPlaneStrain):
  """
  Python object implementing plane strain isotropic linear GenMaxwell
  viscoelastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="genmaxwellplanestrain"):
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
            'data': ["stress_zz_initial",
                     "total_strain", "stress",
                     "viscous_strain_1", 
                     "viscous_strain_2", 
                     "viscous_strain_3",
                     ]}}
    self._loggingPrefix = "MaGM2D "
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleGenMaxwellPlaneStrain.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with GenMaxwellPlaneStrain.
  """
  return GenMaxwellPlaneStrain()


# End of file 
