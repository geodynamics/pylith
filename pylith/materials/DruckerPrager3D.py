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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/DruckerPrager3D.py
##
## @brief Python object implementing 3-D isotropic Drucker-Prager
## elastoplastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import DruckerPrager3D as ModuleDruckerPrager3D

# DruckerPrager3D class
class DruckerPrager3D(ElasticMaterial, ModuleDruckerPrager3D):
  """
  Python object implementing 3-D isotropic Drucker-Prager elastoplastic
  material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="druckerprager3d"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["mu", "lambda", "density", 
                     "alpha_yield", "beta", "alpha_flow"],
            'data': ["total_strain", "stress", "plastic_strain"]}}
    self._loggingPrefix = "MaDP3D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleDruckerPrager3D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with DruckerPrager3D.
  """
  return DruckerPrager3D()


# End of file 
