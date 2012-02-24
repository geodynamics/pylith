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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/ElasticStrain1D.py
##
## @brief Python object implementing 1-D linear elastic material with
## axial strain.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import ElasticStrain1D as ModuleElasticStrain1D

# ElasticStrain1D class
class ElasticStrain1D(ElasticMaterial, ModuleElasticStrain1D):
  """
  Python object implementing 1-D linear elastic material with axial strain.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticstrain1d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["mu", "lambda", "density"],
            'data': ["total_strain", "stress"]}}
    self._loggingPrefix = "MaSt1D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleElasticStrain1D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with ElasticStrain1D.
  """
  return ElasticStrain1D()


# End of file 
