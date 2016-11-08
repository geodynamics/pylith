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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/AuxFieldsIsotropicLinearElasticity.py
##
## @brief Python subfields container for isotropic, linear elasticity
## subfields.

from pylith.utils.PetscComponent import PetscComponent

# AuxFieldsIsotropicLinearElasticity class
class AuxFieldsIsotropicLinearElasticity(PetscComponent):
  """
  Python subfields container for isotropic, linear elasticity subfields.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """Python object for managing AuxFieldsIsotropicLinearElasticity
    facilities and properties.

    """
    
    ## @class Inventory
    ## Python object for managing AuxFieldsIsotropicLinearElasticity facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b density Density subfield.
    ## @li \b bulk_modulus Bulk modulus subfield.
    ## @li \b shear_modulus Shear modulus subfield.

    import pyre.inventory

    from pylith.topology.AuxSubfield import AuxSubfield

    density = pyre.inventory.facility("density", family="auxiliary_subfield", factory=AuxSubfield)
    density.meta['tip'] = "Density subfield."

    shearModulus = pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=AuxSubfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pyre.inventory.facility("bulk_modulus", family="auxiliary_subfield", factory=AuxSubfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    bodyForce = pyre.inventory.facility("body_force", family="auxiliary_subfield", factory=AuxSubfield)
    bodyForce.meta['tip'] = "Body force subfield."

    referenceStress = pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=AuxSubfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=AuxSubfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="auxfieldsisotropiclinearelasticity"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="auxiliary_fields")
    return


  def _configure(self):
    PetscComponent._configure(self)
    self.density = self.inventory.density
    self.shearModulus = self.inventory.shearModulus
    self.bulkModulus = self.inventory.bulkModulus
    self.bodyForce = self.inventory.bodyForce
    self.referenceStress = self.inventory.referenceStress
    self.referenceStrain = self.inventory.referenceStrain
    return


# End of file 
