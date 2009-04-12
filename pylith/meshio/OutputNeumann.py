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

## @file pyre/meshio/OutputNeumann.py
##
## @brief Python object for managing output of finite-element
## information for Neumann boundary conditions.
##
## Factory: output_manager

from OutputManager import SubMeshOutputManager

# OutputNeumann class
class OutputNeumann(SubMeshOutputManager):
  """
  Python object for managing output of finite-element information for
  Neumann boundary conditions.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SubMeshOutputManager.Inventory):
    """
    Python object for managing OutputNeumann facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputNeumann facilities and properties.
    ##
    ## \b Properties
    ## @li \b cell_info_fields Names of cell info fields to output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    cellInfoFields = pyre.inventory.list("cell_info_fields",
                                         default=[])
    cellInfoFields.meta['tip'] = "Names of cell info fields to output."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputneumann"):
    """
    Constructor.
    """
    SubMeshOutputManager.__init__(self, name)
    return

    
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    SubMeshOutputManager._configure(self)
    self.vertexInfoFields = []
    self.vertexDataFields = []
    self.cellInfoFields = self.inventory.cellInfoFields
    self.cellDataFields = []
    return


# FACTORIES ////////////////////////////////////////////////////////////

def submesh_output_manager():
  """
  Factory associated with OutputNeumann.
  """
  return OutputNeumann()


# End of file 
