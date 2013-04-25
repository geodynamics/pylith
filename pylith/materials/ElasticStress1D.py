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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/ElasticStress1D.py
##
## @brief Python object implementing 1-D linear elastic material with
## axial stress.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import ElasticStress1D as ModuleElasticStress1D

# ElasticStress1D class
class ElasticStress1D(ElasticMaterial, ModuleElasticStress1D):
  """
  Python object implementing 1-D linear elastic material with axial stress.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticstress1d"):
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
    self._loggingPrefix = "MaSn1D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleElasticStress1D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticStress1D.
  """
  return ElasticStress1D()


# End of file 
