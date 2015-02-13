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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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
    ## @li \b output Output manager associated with material data.
    ## @li \b db_initial_stress Database for initial stress.
    ## @li \b db_initial_strain Database for initial strain.

    import pyre.inventory

    from pylith.meshio.OutputMatElastic import OutputMatElastic
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputMatElastic)
    output.meta['tip'] = "Output manager for elastic material information."

    from pylith.utils.NullComponent import NullComponent
    dbInitialStress = pyre.inventory.facility("db_initial_stress",
                                              family="spatial_database",
                                              factory=NullComponent)
    dbInitialStress.meta['tip'] = "Database for initial stress."

    dbInitialStrain = pyre.inventory.facility("db_initial_strain",
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

  
  def _modelMemoryUse(self):
    """
    Model allocated memory.
    """
    Material._modelMemoryUse(self)
    self.perfLogger.logFields('Materials', self.initialFields())
    return


# End of file
