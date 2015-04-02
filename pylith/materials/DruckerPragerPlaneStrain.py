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

## @file pylith/materials/DruckerPragerPlaneStrain.py
##
## @brief Python object implementing 2-D plane strain Drucker-Prager
## elastoplastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import DruckerPragerPlaneStrain as ModuleDruckerPragerPlaneStrain

# Validator to fit to Mohr-Coulomb
def validateFitMohrCoulomb(value):
  """
  Validate fit to Mohr-Coulomb yield surface.
  """
  if not value in ["inscribed", "middle", "circumscribed"]:
    raise ValueError("Unknown fit to Mohr-Coulomb yield surface.")
  return value


# DruckerPragerPlaneStrain class
class DruckerPragerPlaneStrain(ElasticMaterial, ModuleDruckerPragerPlaneStrain):
  """
  Python object implementing 2-D plane strain Drucker-Prager elastoplastic
  material.

  Factory: material.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ElasticMaterial.Inventory):
    """
    Python object for managing DruckerPragerPlaneStrain facilities and
    properties.
    """
    
    ## @class Inventory
    ## Python object for managing DruckerPragerPlaneStrain facilities and
    ## properties.
    ##
    ## \b Properties
    ## @li \b fit_mohr_coulomb Fit to Mohr-Coulomb yield surface.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pylith.meshio.OutputMatElastic import OutputMatElastic
    fitMohrCoulomb = pyre.inventory.str("fit_mohr_coulomb", default="inscribed",
                                        validator=validateFitMohrCoulomb)
    fitMohrCoulomb.meta['tip'] = "Fit to Mohr-Coulomb yield surface."

    allowTensileYield = pyre.inventory.bool("allow_tensile_yield", default=False)
    allowTensileYield.meta['tip'] = "Extend yield surface past tip of cone to allow yielding with tensile stresses."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="druckerpragerplanestrain"):
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
                     "alpha_yield", "beta", "alpha_flow"],
            'data': ["total_strain", "stress", "stress_zz_initial",
                     "plastic_strain"]}}
    self._loggingPrefix = "MaDP2D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    ElasticMaterial._configure(self)
    if self.inventory.fitMohrCoulomb == "inscribed":
      fitEnum = ModuleDruckerPragerPlaneStrain.MOHR_COULOMB_INSCRIBED
    elif self.inventory.fitMohrCoulomb == "middle":
      fitEnum = ModuleDruckerPragerPlaneStrain.MOHR_COULOMB_MIDDLE
    elif self.inventory.fitMohrCoulomb == "circumscribed":
      fitEnum = ModuleDruckerPragerPlaneStrain.MOHR_COULOMB_CIRCUMSCRIBED
    else:
      raise ValueError("Unknown fit to Mohr-Coulomb yield surface.")
    ModuleDruckerPragerPlaneStrain.fitMohrCoulomb(self, fitEnum)
    ModuleDruckerPragerPlaneStrain.allowTensileYield(self, self.inventory.allowTensileYield)
    return

  
  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleDruckerPragerPlaneStrain.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with DruckerPragerPlaneStrain.
  """
  return DruckerPragerPlaneStrain()


# End of file 
