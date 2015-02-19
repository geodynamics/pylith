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

## @file pylith/materials/PowerLawPlaneStrain.py
##
## @brief Python object implementing 2-D plane strain power-law
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import PowerLawPlaneStrain as ModulePowerLawPlaneStrain

# PowerLawPlaneStrain class
class PowerLawPlaneStrain(ElasticMaterial, ModulePowerLawPlaneStrain):
  """
  Python object implementing 2-D plane strain power-law viscoelastic material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="powerlawplanestrain"):
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
            'data': ["total_strain", "stress",
                     "stress_zz_initial", "stress4", "viscous_strain"]}}
    self._loggingPrefix = "MaPL2D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModulePowerLawPlaneStrain.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with PowerLawPlaneStrain.
  """
  return PowerLawPlaneStrain()


# End of file 
