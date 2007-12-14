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

## @file pyre/meshio/OutputManager.py
##
## @brief Python abstract base class for managing output of
## finite-element information.
##
## Factory: output_manager

from pyre.components.Component import Component

# OutputManager class
class OutputManager(Component):
  """
  Python abstract base class for managing output of finite-element
  information.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing OutputManager facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing OutputManager facilities and properties.
    ##
    ## \b Properties
    ## @li \b vertex_fields Fields of vertex data to output.
    ## @li \b cell_fields Fields of cell data to output.
    ## @li \b output_freq Flag indicating whether to use 'time_step' or 'skip'
    ##   to set frequency of solution output.
    ## @li \b time_step Time step between solution output.
    ## @li \b skip Number of time steps to skip between solution output.
    ##
    ## \b Facilities
    ## @li \b writer Writer for data.
    ## @li \b vertex_filter Filter for vertex data.
    ## @li \b cell_filter Filter for cell data.

    import pyre.inventory

    vertexFields = pyre.inventory.list("vertex_fields", default=[])
    vertexFields.meta['tip'] = "Fields of vertex data to output."

    cellFields = pyre.inventory.list("cell_fields", default=[])
    cellFields.meta['tip'] = "Fields of cell data to output."

    writer = pyre.inventory.facility("writer", factory=DataWriterVTK,
                                     family="data_writer")
    writer.meta['tip'] = "Writer for data."

    vertexFilter = pyre.inventory.facility("vertex_filter",
                                           factory=VertexFilter,
                                           family="output_vertex_filter")
    vertexFilter.meta['tip'] = "Filter for vertex data."
                                     
    cellFilter = pyre.inventory.facility("cell_filter",
                                           factory=CellFilter,
                                           family="output_cell_filter")
    cellFilter.meta['tip'] = "Filter for cell data."
                                     

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputmanager"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="outputmanager")
    self._loggingPrefix = "OutM "
    self._t = None
    self._istep = None
    self.cppHandle = None
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self, quadrature):
    """
    Initialize output manager.
    """
    self.cellFilter.initialize(quadrature)
    return


  def open(self, mesh):
    """
    Prepare for writing data.
    """
    self._setupLogging()
    logEvent = "%sopen" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    self.mesh = mesh # Keep handle to mesh
    self._sync()
    
    assert(None != self.cppHandle)
    assert(None != mesh.cppHandle)
    assert(None != mesh.coordsys.cppHandle)
    self.cppHandle.open(mesh.cppHandle, mesh.coordsys.cppHandle)

    self._logger.eventEnd(logEvent)    
    return


  def close(self):
    """
    Perform post-write cleanup.
    """
    logEvent = "%sclose" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    assert(None != self.cppHandle)
    self.cppHandle.close()

    self._logger.eventEnd(logEvent)    
    return


  def writeFields(self, t, istep, fields=None):
    """
    Write vertex and cell fields.

    @param fields FieldsManager containing fields (if not in mesh).
    """
    logEvent = "%swriteVertex" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    write = False
    if self.istep == None or not "value" in dir(self._t):
      write = True
    elif self.outputFreq == "skip":
      if istep > self._istep + self.skip:
        write = True
    elif t >= self._t + self._dt:
      write = True
    if write:
      self._info.log("Writing fields.")
      assert(None != self.cppHandle)
      assert(None != mesh.cppHandle)
      assert(None != mesh.coordsys.cppHandle)
      self.cppHandle.writeFields(t, fields,
                                 self.mesh.cppHandle,
                                 self.mesh.coordsys.cppHandle)
      self._istep = istep
      self._t = t

    self._logger.eventEnd(logEvent)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.vertexFields = self.inventory.vertexFields
    self.cellFields = self.inventory.cellFields
    self.writer = self.inventory.writer
    self.vertexFilter = self.inventory.vertexFilter
    self.cellFilter = self.inventory.cellFilter
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    if None == self.cppHandle:
      import pylith.meshio.meshio as bindings
      self.cppHandle = bindings.OutputManager()
    self.cppHandle.vertexFields = self.vertexFields
    self.cppHandle.cellFields = self.cellFields
    self.cppHandle.vertexFilter = self.vertexFilter.cppHandle
    self.cppHandle.cellFilter = self.cellFilter.cppHandle
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("FE Output")
    logger.initialize()

    events = ["open",
              "close",
              "writeFields"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputManager()


# End of file 
