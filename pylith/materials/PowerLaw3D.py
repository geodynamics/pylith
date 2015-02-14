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

## @file pylith/materials/PowerLaw3D.py
##
## @brief Python object implementing 3-D isotropic power-law
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import PowerLaw3D as ModulePowerLaw3D

# PowerLaw3D class
class PowerLaw3D(ElasticMaterial, ModulePowerLaw3D):
  """
  Python object implementing 3-D isotropic power-law viscoelastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="powerlaw3d"):
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
                     "reference_strain_rate", "reference_stress",
                     "power_law_exponent"],
            'data': ["total_strain", "stress", "viscous_strain"]}}
    self._loggingPrefix = "MaPL3D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModulePowerLaw3D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with PowerLaw3D.
  """
  return PowerLaw3D()


# End of file 
