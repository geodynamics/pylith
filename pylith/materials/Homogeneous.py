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

## @file pylith/materials/Homogeneous.py
## @brief Python materials container with one material.

from MaterialsBin import MaterialsBin

# Homogeneous class
class Homogeneous(MaterialsBin):
  """Python materials container with one material."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MaterialsBin.Inventory):
    """Python object for managing Homogeneous facilities and properties."""
    
    ## @class Inventory
    ## Python object for managing Homogeneous facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b material Material in problem

    import pyre.inventory

    from ElasticIsotropic import ElasticIsotropic
    material = pyre.inventory.facility("material", factory=ElasticIsotropic)
    material.meta['tip'] = "Material in problem."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="Homogeneous"):
    """Constructor."""
    
    MaterialsBin.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set attributes from inventory."""

    MaterialsBin._configure(self)
    self.materials = [self.inventory.material]
    return

  
 # version
__id__ = "$Id$"

# End of file 
