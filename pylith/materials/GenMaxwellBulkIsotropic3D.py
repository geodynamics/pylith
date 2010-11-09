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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/GenMaxwellBulkIsotropic3D.py
##
## @brief Python object implementing 3-D isotropic linear GenMaxwellBulk
## viscoelastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import GenMaxwellBulkIsotropic3D as ModuleGenMaxwellBulkIsotropic3D

# GenMaxwellBulkIsotropic3D class
class GenMaxwellBulkIsotropic3D(ElasticMaterial, ModuleGenMaxwellBulkIsotropic3D):
  """
  Python object implementing 3-D isotropic linear GenMaxwellBulk viscoelastic
  material.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="genmaxwellbulkisotropic3d"):
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
                     "shear_ratio", "bulk_ratio"
                     "maxwell_time_shear",
                     "maxwell_time_bulk"],
            'data': ["total_strain", "stress",
                     "viscous_strain_1", 
                     "viscous_strain_2", 
                     "viscous_strain_3",
                     "viscous_strain_1_bulk", 
                     "viscous_strain_2_bulk", 
                     "viscous_strain_3_bulk",
                     ]}}
    self._loggingPrefix = "MaMx3D "
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleGenMaxwellBulkIsotropic3D.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with GenMaxwellBulkIsotropic3D.
  """
  return GenMaxwellBulkIsotropic3D()


# End of file 
