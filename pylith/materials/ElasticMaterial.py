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
    ## @li \b initial_stress_db Database for initial stress.
    ## @li \b initial_strain_db Database for initial strain.

    import pyre.inventory

    from pylith.meshio.OutputMatElastic import OutputMatElastic
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputMatElastic)
    output.meta['tip'] = "Output manager for elastic material information."

    from pylith.utils.NullComponent import NullComponent
    dbInitialStress = pyre.inventory.facility("initial_stress_db",
                                              family="spatial_database",
                                              factory=NullComponent)
    dbInitialStress.meta['tip'] = "Database for initial stress."

    dbInitialStrain = pyre.inventory.facility("initial_strain_db",
                                              family="spatial_database",
                                              factory=NullComponent)
    dbInitialStrain.meta['tip'] = "Database for initial strain."


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
    from pylith.utils.NullComponent import NullComponent
    if not isinstance(self.inventory.dbInitialStress, NullComponent):
      self.dbInitialStress(self.inventory.dbInitialStress)
    if not isinstance(self.inventory.dbInitialStrain, NullComponent):
      self.dbInitialStrain(self.inventory.dbInitialStrain)
    return

  
# End of file 
