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

## @file pylith/materials/MaxwellPlaneStrain.py
##
## @brief Python object implementing plane strain linear Maxwell
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import MaxwellPlaneStrain as ModuleMaxwellPlaneStrain

# MaxwellPlaneStrain class
class MaxwellPlaneStrain(ElasticMaterial, ModuleMaxwellPlaneStrain):
  """
  Python object implementing plane strain linear Maxwell viscoelastic
  material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="maxwellplanestrain"):
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
            'data': ["total_strain", "stress",
                     "stress_zz_initial", "viscous_strain"]}}
    self._loggingPrefix = "MaMx2D "
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleMaxwellPlaneStrain.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with MaxwellPlaneStrain.
  """
  return MaxwellPlaneStrain()


# End of file 
