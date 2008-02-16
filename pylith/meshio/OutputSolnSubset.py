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

## @file pyre/meshio/OutputSolnSubset.py
##
## @brief Python object for managing output of finite-element solution
## information over a subdomain.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputSolnSubset class
class OutputSolnSubset(OutputManager):
  """
  Python object for managing output of finite-element solution
  information over a subdomain.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(OutputManager.Inventory):
    """
    Python object for managing OutputSolnSubset facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputSolnSubset facilities and properties.
    ##
    ## \b Properties
    ## @li \b vertex_data_fields Names of vertex data fields to output.
    ## @li \b label Name identifier for subdomain.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                           default=["displacements"])
    vertexDataFields.meta['tip'] = "Names of vertex data fields to output."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Label identifier for subdomain."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputsolnsubset"):
    """
    Constructor.
    """
    OutputManager.__init__(self, name)
    self.cppHandle = None
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': ["displacements"]},
         'cell': \
           {'info': [],
            'data': []}}
    return


  def preinitialize(self):
    """
    Do
    """
    OutputManager.preinitialize(self, dataProvider=self)
    import meshio as bindings
    self.cppHandle = bindings.OutputSolnSubset()
    self.cppHandle.label = self.label
    return
  

  def initialize(self, mesh):
    """
    Initialize output manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    from pylith.topology.Mesh import Mesh
    self.mesh = Mesh()
    self.mesh.initialize(mesh.coordsys)
    assert(None != self.cppHandle)
    self.cppHandle.mesh(self.mesh.cppHandle, mesh.cppHandle)

    self.mesh.view()

    OutputManager.initialize(self)

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


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManager._configure(self)
    self.vertexDataFields = self.inventory.vertexDataFields
    self.label = self.inventory.label
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputSolnSubset.
  """
  return OutputSolnSubset()


# End of file 
