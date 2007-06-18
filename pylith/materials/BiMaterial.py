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

## @file pylith/materials/BiMaterial.py
##
## @brief Python materials container with one material.
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# BiMaterial class
class BiMaterial(ObjectBin):
  """
  Python materials container with two materials.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
    """
    Python object for managing BiMaterial facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BiMaterial facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b one One material in problem
    ## @li \b two Other material in problem.

    import pyre.inventory

    from ElasticIsotropic3D import ElasticIsotropic3D

    one = pyre.inventory.facility("one", family="material",
                                  factory=ElasticIsotropic3D)
    one.meta['tip'] = "One material in problem."
    
    two = pyre.inventory.facility("two", family="material",
                                  factory=ElasticIsotropic3D)
    two.meta['tip'] = "Other material in problem."


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
    self.bin = [self.inventory.one,
                self.inventory.two]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with BiMaterial.
  """
  return BiMaterial()


# End of file 
