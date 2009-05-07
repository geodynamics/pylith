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
    ## @li \b use_initial_stress Use initial stress (true) or not (false).
    ## @li \b use_initial_strain Use initial strain (true) or not (false).
    ##
    ## \b Facilities
    ## @li \b output Output manager associated with fault data.
    ## @li \b initial_stress_db Database for initial stress.
    ## @li \b initial_strain_db Database for initial strain.

    import pyre.inventory

    useInitialStress = pyre.inventory.bool("use_initial_stress", default=False)
    useInitialStress.meta['tip'] = "Use initial stress for material."

    useInitialStrain = pyre.inventory.bool("use_initial_strain", default=False)
    useInitialStrain.meta['tip'] = "Use initial strain for material."

    from pylith.meshio.OutputMatElastic import OutputMatElastic
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputMatElastic)
    output.meta['tip'] = "Output manager for elastic material information."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbInitialStress = pyre.inventory.facility("initial_stress_db",
                                              family="spatial_database",
                                              factory=SimpleDB)
    dbInitialStress.meta['tip'] = "Database for initial stress."

    dbInitialStrain = pyre.inventory.facility("initial_strain_db",
                                              family="spatial_database",
                                              factory=SimpleDB)
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
    if self.inventory.useInitialStress:
      self.dbInitialStress(self.inventory.dbInitialStress)
    if self.inventory.useInitialStrain:
      self.dbInitialStrain(self.inventory.dbInitialStrain)
    return

  
# End of file 
