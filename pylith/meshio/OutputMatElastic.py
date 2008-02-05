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

## @file pyre/meshio/OutputMatElastic.py
##
## @brief Python object for managing output of finite-element
## information for faults with kinematic ruptures.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputMatElastic class
class OutputMatElastic(OutputManager):
  """
  Python object for managing output of finite-element information for
  faults with kinematic ruptures.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(OutputManager.Inventory):
    """
    Python object for managing OutputMatElastic facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputMatElastic facilities and properties.
    ##
    ## \b Properties
    ## @li \b vertex_info_fields Names of vertex info fields to output.
    ## @li \b vertex_data_fields Names of vertex data fields to output.
    ## @li \b cell_info_fields Names of cell info fields to output.
    ## @li \b cell_data_fields Names of cell data fields to output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    vertexInfoFields = pyre.inventory.list("vertex_info_fields",
                                           default=[])
    vertexInfoFields.meta['tip'] = "Names of vertex info fields to output."

    vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                           default=[])
    vertexDataFields.meta['tip'] = "Names of vertex data fields to output."

    cellInfoFields = pyre.inventory.list("cell_info_fields", default=[])
    cellInfoFields.meta['tip'] = "Names of cell info fields to output."

    cellDataFields = pyre.inventory.list("cell_data_fields", 
                                         default=[])
    cellDataFields.meta['tip'] = "Names of cell data fields to output."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputmanager"):
    """
    Constructor.
    """
    OutputManager.__init__(self, name)
    return

    
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManager._configure(self)
    self.vertexInfoFields = self.inventory.vertexInfoFields
    self.vertexDataFields = self.inventory.vertexDataFields
    self.cellInfoFields = self.inventory.cellInfoFields
    self.cellDataFields = self.inventory.cellDataFields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputMatElastic.
  """
  return OutputMatElastic()


# End of file 
