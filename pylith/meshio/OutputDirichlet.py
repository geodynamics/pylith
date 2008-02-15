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

## @file pyre/meshio/OutputDirichlet.py
##
## @brief Python object for managing output of finite-element
## information for Dirichlet boundary conditions.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputDirichlet class
class OutputDirichlet(OutputManager):
  """
  Python object for managing output of finite-element information for
  Dirichlet boundary conditions.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(OutputManager.Inventory):
    """
    Python object for managing OutputDirichlet facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputDirichlet facilities and properties.
    ##
    ## \b Properties
    ## @li \b vertex_info_fields Names of vertex info fields to output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    vertexInfoFields = pyre.inventory.list("vertex_info_fields",
                                           default=[])
    vertexInfoFields.meta['tip'] = "Names of vertex info fields to output."
    

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
    self.vertexDataFields = []
    self.cellInfoFields = []
    self.cellDataFields = []
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputDirichlet.
  """
  return OutputDirichlet()


# End of file 
