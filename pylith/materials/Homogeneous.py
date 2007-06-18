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
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# Homogeneous class
class Homogeneous(ObjectBin):
  """
  Python materials container with one material.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
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
    ObjectBin.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set attributes from inventory.
    """
    ObjectBin._configure(self)
    self.bin = [self.inventory.material]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with Homogeneous.
  """
  return Homogeneous()


# End of file 
