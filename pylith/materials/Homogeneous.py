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
##
## @brief Python materials container with one material.

from pyre.components.Component import Component

# Homogeneous class
class Homogeneous(Component):
  """
  Python materials container with one material.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Homogeneous facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Homogeneous facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b material Material in problem

    import pyre.inventory

    from ElasticIsotropic3D import ElasticIsotropic3D
    material = pyre.inventory.facility("material", family="material",
                                       factory=ElasticIsotropic3D)
    material.meta['tip'] = "Material in problem."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="homogeneous"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return


# End of file 
