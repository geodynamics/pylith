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

## @file pyre/meshio/OutputSoln.py
##
## @brief Python object for managing output of finite-element
## solution information.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputSoln class
class OutputSoln(OutputManager):
  """
  Python object for managing output of finite-element solution
  information.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(OutputManager.Inventory):
    """
    Python object for managing OutputSoln facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputSoln facilities and properties.
    ##
    ## \b Properties
    ## @li \b vertex_data_fields Names of vertex data fields to output.
    ## @li \b cell_info_fields Names of cell info fields to output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                           default=["displacements"])
    vertexDataFields.meta['tip'] = "Names of vertex data fields to output."

    cellInfoFields = pyre.inventory.list("cell_info_fields", default=[])
    cellInfoFields.meta['tip'] = "Names of cell info fields to output."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputsoln"):
    """
    Constructor.
    """
    OutputManager.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': ["displacements"]},
         'cell': \
           {'info': ["replaced_cells"],
            'data': []}}
    return


  def preinitialize(self):
    """
    Do
    """
    OutputManager.preinitialize(self, dataProvider=self)
    return
  

  def initialize(self, mesh, normalizer):
    """
    Initialize output manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    self.mesh = mesh
    OutputManager.initialize(self, normalizer)

    self._logger.eventEnd(logEvent)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.mesh, None, None)


  def getVertexField(self, name, fields):
    """
    Get vertex field.
    """
    field = None
    fieldType = None
    if name == "displacements":
      field = fields.getSolution()
      fieldType = 1 # vector field
    else:
      raise ValueError, "Vertex field '%s' not available." % name
    return (field, fieldType)


  def getCellField(self, name):
    """
    Get vertex field.
    """
    field = None
    fieldType = None
    if name == "replaced_cells":
      field = self.mesh.getRealSection("replaced_cells")
      fieldType = 0 # scalar field
    else:
      raise ValueError, "Vertex field '%s' not available." % name
    return (field, fieldType)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManager._configure(self)
    self.vertexDataFields = self.inventory.vertexDataFields
    self.cellInfoFields = self.inventory.cellInfoFields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputSoln.
  """
  return OutputSoln()


# End of file 
