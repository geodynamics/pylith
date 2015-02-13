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

## @file pylith/materials/DruckerPrager3D.py
##
## @brief Python object implementing 3-D isotropic Drucker-Prager
## elastoplastic material.
##
## Factory: material.

from ElasticMaterial import ElasticMaterial
from materials import DruckerPrager3D as ModuleDruckerPrager3D

# Validator to fit to Mohr-Coulomb
def validateFitMohrCoulomb(value):
  """
  Validate fit to Mohr-Coulomb yield surface.
  """
  if not value in ["inscribed", "middle", "circumscribed"]:
    raise ValueError("Unknown fit to Mohr-Coulomb yield surface.")
  return value


# DruckerPrager3D class
class DruckerPrager3D(ElasticMaterial, ModuleDruckerPrager3D):
  """
  Python object implementing 3-D isotropic Drucker-Prager elastoplastic
  material.

  Factory: material.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ElasticMaterial.Inventory):
    """
    Python object for managing FaultCohesiveKin facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FaultCohesiveKin facilities and properties.
    ##
    ## \b Properties
    ## @li \b fit_mohr_coulomb Fit to Mohr-Coulomb yield surface.
    ## @li \b allow_tensile_yield If true, allow yield beyond tensile strength;
    ##   otherwise an exception occurs for excessive tensile sttress.
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
           {'info': ["mu", "lambda", "density", "stable_dt_implicit", "stable_dt_explicit",
                     "alpha_yield", "beta", "alpha_flow"],
            'data': ["total_strain", "stress", "plastic_strain"]}}
    self._loggingPrefix = "MaDP3D "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    ElasticMaterial._configure(self)
    if self.inventory.fitMohrCoulomb == "inscribed":
      fitEnum = ModuleDruckerPrager3D.MOHR_COULOMB_INSCRIBED
    elif self.inventory.fitMohrCoulomb == "middle":
      fitEnum = ModuleDruckerPrager3D.MOHR_COULOMB_MIDDLE
    elif self.inventory.fitMohrCoulomb == "circumscribed":
      fitEnum = ModuleDruckerPrager3D.MOHR_COULOMB_CIRCUMSCRIBED
    else:
      raise ValueError("Unknown fit to Mohr-Coulomb yield surface.")
    ModuleDruckerPrager3D.fitMohrCoulomb(self, fitEnum)
    ModuleDruckerPrager3D.allowTensileYield(self, self.inventory.allowTensileYield)
    return

  
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
