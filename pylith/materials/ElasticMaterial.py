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

## @file pylith/materials/ElasticMaterial.py
##
## @brief Python abstract base class for managing physical properties
## of an elastic material.
##
## Factory: material

from Material import Material

# ElasticMaterial class
class ElasticMaterial(Material):
  """
  Python abstract base class for managing physical properties of an
  elastic material.

  Factory: material
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Material.Inventory):
    """
    Python object for managing FaultCohesiveKin facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FaultCohesiveKin facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b output Output manager associated with fault data.

    import pyre.inventory

    from pylith.meshio.OutputMatElastic import OutputMatElastic
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputMatElastic)
    output.meta['tip'] = "Output manager for elastic material information."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticmaterial"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Material._configure(self)
    self.output = self.inventory.output
    return

  
# End of file 
