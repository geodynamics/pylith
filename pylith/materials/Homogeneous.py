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
## Factory: materials_bin

from MaterialsBin import MaterialsBin

# Homogeneous class
class Homogeneous(MaterialsBin):
  """
  Python materials container with one material.

  Factory: materials_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MaterialsBin.Inventory):
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
    MaterialsBin.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set attributes from inventory.
    """
    MaterialsBin._configure(self)
    self.materials = [self.inventory.material]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def materials_bin():
  """
  Factory associated with Homogeneous.
  """
  return Homogeneous()


# End of file 
