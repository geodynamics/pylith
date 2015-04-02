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

## @file pylith/materials/ElasticPlaneStrain.py
##
## @brief Python object implementing 1-D isotropic linear elastic
## material for plane strain.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import ElasticPlaneStrain as ModuleElasticPlaneStrain

# ElasticPlaneStrain class
class ElasticPlaneStrain(ElasticMaterial, ModuleElasticPlaneStrain):
  """
  Python object implementing 2-D isotropic linear elastic material for
  plane strain.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticplanestrain"):
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
    self._loggingPrefix = "MaPlSn "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleElasticPlaneStrain.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticPlaneStrain.
  """
  return ElasticPlaneStrain()


# End of file 
