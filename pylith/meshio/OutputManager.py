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
## @brief Python object for managing output of finite-element
## information.
##
## Factory: output_manager

from pyre.components.Component import Component

# OutputManager class
class OutputManager(Component):
  """
  Python object for managing output of finite-element information.

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
    ## @li \b vertex_fields Names of vertex fields to output.
    ## @li \b cell_fields Names of cell fields to output.
    ## @li \b output_freq Flag indicating whether to use 'time_step' or 'skip'
    ##   to set frequency of solution output.
    ## @li \b time_step Time step between solution output.
    ## @li \b skip Number of time steps to skip between solution output.
    ##
    ## \b Facilities
    ## @li \b writer Writer for data.
    ## @li \b coordsys Coordinate system for output.
    ## @li \b vertex_filter Filter for vertex data.
    ## @li \b cell_filter Filter for cell data.

    import pyre.inventory

    outputFreq = pyre.inventory.str("output_freq", default="skip",
             validator=pyre.inventory.choice(["skip", "time_step"]))
    outputFreq.meta['tip'] = "Flag indicating whether to use 'time_step' " \
                             "or 'skip' to set frequency of output."

    from pyre.units.time import s
    dt = pyre.inventory.dimensional("time_step", default=1.0*s)
    dt.meta['tip'] = "Time step between output."

    skip = pyre.inventory.int("skip", default=0,
                              validator=pyre.inventory.greaterEqual(0))
    skip.meta['tip'] = "Number of time steps to skip between output."

    vertexFields = pyre.inventory.list("vertex_fields", default=[])
    vertexFields.meta['tip'] = "Fields of vertex data to output."

    cellFields = pyre.inventory.list("cell_fields", default=[])
    cellFields.meta['tip'] = "Fields of cell data to output."

    from DataWriterVTK import DataWriterVTK
    writer = pyre.inventory.facility("writer", factory=DataWriterVTK,
                                     family="data_writer")
    writer.meta['tip'] = "Writer for data."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system for output."
  
    from VertexFilter import VertexFilter
    vertexFilter = pyre.inventory.facility("vertex_filter",
                                           factory=VertexFilter,
                                           family="output_vertex_filter")
    vertexFilter.meta['tip'] = "Filter for vertex data."
                                     
    from CellFilter import CellFilter
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
    self.cppHandle = None
    self.coordsys = None
    self.mesh = None
    self._t = None
    self._istep = None
    self.vertexFields = None
    self.cellFields = None
    self._fieldTranslator = copyTranslator
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    # Verify fields requested for output are available by creating map
    # of names of requested fields to mesh labels.
    self._createFieldDicts()
    return


  def fieldTranslator(self, translator):
    """
    Set function to call to translate names of fields to mesh labels.
    """
    self._fieldTranslator = translator
    return


  def initialize(self, quadrature=None):
    """
    Initialize output manager.
    """
    self._createFieldDicts()
    
    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for output is unknown."
    self.coordsys.initialize()

    self.cellFilter.initialize(quadrature)
    self.writer.initialize()
    self._sync()
    return


  def open(self, mesh):
    """
    Prepare for writing data.
    """
    self._setupLogging()
    logEvent = "%sopen" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    self.mesh = mesh # Keep handle to mesh
    
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


  def openTimeStep(self, t, istep):
    """
    Prepare for writing solution to file.
    """
    logEvent = "%sopenStep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    
    self._info.log("Preparing for writing solution to file.")

    write = False
    if self._istep == None or not "value" in dir(self._t):
      write = True
    elif self.outputFreq == "skip":
      if istep > self._istep + self.skip:
        write = True
    elif t >= self._t + self.dt:
      write = True
    if write:
      self._istep = istep
      self._t = t
    self.writeFlag = write

    assert(self.cppHandle != None)
    assert(self.mesh.cppHandle != None)
    assert(self.mesh.coordsys.cppHandle != None)
    self.cppHandle.openTimeStep(t.value,
                                self.mesh.cppHandle,
                                self.mesh.coordsys.cppHandle)

    self._logger.eventEnd(logEvent)    
    return


  def closeTimeStep(self):
    """
    Cleanup after writing solution to file.
    """
    logEvent = "%scloseStep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    
    self._info.log("Cleaning up afterwriting solution to file.")

    self.writeFlag = False

    assert(self.cppHandle != None)
    self.cppHandle.closeTimeStep()

    self._logger.eventEnd(logEvent)    
    return


  def appendVertexField(self, t, istep, name, field, dim=0):
    """
    Write field over vertices at time t to file.
    """
    logEvent = "%swriteVertex" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    if self.writeFlag:
      self._info.log("Writing solution field '%s'." % name)
      assert(self.cppHandle != None)
      assert(self.mesh.cppHandle != None)
      self.cppHandle.appendVertexField(t.value, name, field,
                                       self.mesh.cppHandle, dim)
      self.istep = istep
      self.t = t

    self._logger.eventEnd(logEvent)
    return


  def appendCellField(self, t, istep, name, field, dim=0):
    """
    Write field over cells at time t to file.
    """
    logEvent = "%swriteCell" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    if self.writeFlag:
      self._info.log("Writing solution field '%s'." % name)
      assert(self.cppHandle != None)
      assert(self.mesh.cppHandle != None)
      self.cppHandle.appendCellField(t.value, name, field, 
                                     self.mesh.cppHandle, dim)
      self.istep = istep
      self.t = t

    self._logger.eventEnd(logEvent)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.outputFreq = self.inventory.outputFreq
    self.dt = self.inventory.dt
    self.skip = self.inventory.skip
    self.coordsys = self.inventory.coordsys
    self.writer = self.inventory.writer
    self.vertexFieldNames = self.inventory.vertexFields
    self.cellFieldNames = self.inventory.cellFields
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

    assert(self.coordsys.cppHandle != None)
    assert(self.writer.cppHandle != None)
      
    self.cppHandle.coordsys = self.coordsys.cppHandle
    self.cppHandle.writer = self.writer.cppHandle
    self.cppHandle.vertexFields = self.vertexFields
    self.cppHandle.cellFields = self.cellFields
    self.cppHandle.vertexFilter = self.vertexFilter.cppHandle
    self.cppHandle.cellFilter = self.cellFilter.cppHandle
    return


  def _createFieldDicts(self):
    """
    Create dictionaries with field names and mesh labels of fields.
    """
    if None != self.vertexFields and None != self.cellFields:
      return
    
    self.vertexFields = {}
    self.cellFields = {}
    for name in self.vertexFieldNames:
      self.vertexFields[name] = self._fieldTranslator(name)
    for name in self.cellFieldNames:
      self.cellFields[name] = self._fieldTranslator(name)
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
              "writeFields",
              "openStep",
              "closeStep",
              "writeVertex",
              "writeCell"]
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


# MISCELLANEOUS ////////////////////////////////////////////////////////

def copyTranslator(name):
  return name


# End of file 
