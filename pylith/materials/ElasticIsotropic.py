#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/ElasticIsotropic.py
## @brief Python object for isotropic linear elastic constitutive model.

from MaterialModel import MaterialModel

# ElasticIsotropic class
class ElasticIsotropic(MaterialModel):
  """Python object for isotropic linear elastic constitutive model."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MaterialModel.Inventory):
    """Python object for managing ElasticIsotropic facilities and
    properties."""

    ## @class Inventory
    ## Python object for managing ElasticIsotropic facilities and properties.
    ##
    ## \b Properties
    ## @li \b mu Lame's constant (shear modulus)
    ## @li \b lambda Lame's constant
    ## @li \b density Mass density
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pyre.units.pressure import GPa
    muLame = pyre.inventory.dimensional("mu", default=3.0*GPa)
    muLame.meta['tip'] = "Lame's constant (shear modulus)."

    lambdaLame = pyre.inventory.dimensional("lambda", default=3.0*GPa)
    lambdaLame.meta['tip'] = "Lame's constant."

    from pyre.units.mass import kg
    from pyre.units.length import m
    density = pyre.inventory.dimensional("density", default=3000.0*kg/m**3)
    density.meta['tip'] = "Mass density."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def elasticityConsts(self):
    size = 36
    values = size*[0]
    return values

  def __init__(self, name="matsolidlinelast"):
    """Constructor."""
    MaterialModel.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set members using inventory."""
    self.muLame = self.inventory.muLame
    self.lambdaLame = self.inventory.lambdaLame
    self.density = self.inventory.density
    return


 # version
__id__ = "$Id$"

# End of file 
