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

from pylith.utils.PetscComponent import PetscComponent
from meshio import MeshOutputManager as ModuleMeshObject
from meshio import SubMeshOutputManager as ModuleSubMeshObject

# OutputManager class
class OutputManager(PetscComponent):
  """
  Python abstract base class for managing output of finite-element
  information.

  \b Properties
  @li \b output_freq Flag indicating whether to use 'time_step' or 'skip'
  to set frequency of solution output.
  @li \b time_step Time step between solution output.
  @li \b skip Number of time steps to skip between solution output.
  
  \b Facilities
  @li \b coordsys Coordinate system for output.
  @li \b vertex_filter Filter for vertex data.
  @li \b cell_filter Filter for cell data.
  """

  # INVENTORY //////////////////////////////////////////////////////////

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
    PetscComponent.__init__(self, name, facility="outputmanager")
    self._loggingPrefix = "OutM "
    self._stepCur = 0
    self._stepWrite = None
    self._tWrite = None
    self.dataProvider = None
    self.vertexInfoFields = []
    self.vertexDataFields = []
    self.cellInfoFields = []
    self.cellDataFields = []

    self._createModuleObj()
    return


  def preinitialize(self, dataProvider):
    """
    Setup output manager.
    """
    self._setupLogging()
    self.dataProvider = dataProvider
    return
  

  def verifyConfiguration(self, mesh):
    """
    Verify compatibility of configuration.
    """
    if None == self.dataProvider:
      raise ValueError("Need to set 'dataProvider' in OutputManager.")
    self._verifyFields(self.dataProvider.availableFields)

    if not "getDataMesh" in dir(self.dataProvider):
      raise TypeError("Data provider must have a 'getDataMesh' function.")

    if len(self.vertexInfoFields) > 0 or len(self.vertexDataFields) > 0:
      if not "getVertexField" in dir(self.dataProvider):
        raise TypeError("Data provider must have a 'getVertexField' function.")

    if len(self.cellInfoFields) > 0 or len(self.cellDataFields) > 0:
      if not "getCellField" in dir(self.dataProvider):
        raise TypeError("Data provider must have a 'getCellField' function.")
    return


  def initialize(self, normalizer, quadrature=None):
    """
    Initialize output manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    # Nondimensionalize time step
    lengthScale = normalizer.timeScale()
    self.dt = normalizer.nondimensionalize(self.dt, lengthScale)

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for output is unknown."
    self.coordsys.initialize()

    self.cellFilter.initialize(quadrature)
    self.writer.initialize(normalizer)

    self._logger.eventEnd(logEvent)
    return


  def open(self, totalTime, numTimeSteps):
    """
    Prepare for writing data.
    """
    logEvent = "%sopen" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    nsteps = numTimeSteps
    if numTimeSteps > 0 and self.outputFreq == "skip" and self.skip > 0:
      nsteps = int(numTimeSteps / (1+self.skip))
    elif numTimeSteps > 0 and self.outputFreq == "time_step":
      nsteps = 1 + int(totalTime / self.dt)

    (mesh, label, labelId) = self.dataProvider.getDataMesh()
    self._open(mesh, nsteps, label, labelId)

    self._logger.eventEnd(logEvent)    
    return


  def close(self):
    """
    Perform post-write cleanup.
    """
    logEvent = "%sclose" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    self._close()

    self._logger.eventEnd(logEvent)    
    return


  def writeInfo(self):
    """
    Write information fields.
    """
    logEvent = "%swriteInfo" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    if len(self.vertexInfoFields) > 0 or len(self.cellInfoFields) > 0:
      t = 0.0
      self.open(totalTime=0.0, numTimeSteps=0)
      (mesh, label, labelId) = self.dataProvider.getDataMesh()
      self._openTimeStep(t, mesh, label, labelId)

      for name in self.vertexInfoFields:
        field = self.dataProvider.getVertexField(name)
        self._appendVertexField(t, field)

      for name in self.cellInfoFields:
        field = self.dataProvider.getCellField(name)
        self._appendCellField(t, field, label, labelId)

      self._closeTimeStep()
      self._close()

    self._logger.eventEnd(logEvent)
    return


  def writeData(self, t, fields):
    """
    Write fields at current time step.
    """
    logEvent = "%swriteData" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    if self._checkWrite(t) and \
           ( len(self.vertexDataFields) > 0 or \
             len(self.cellDataFields) ) > 0:

      (mesh, label, labelId) = self.dataProvider.getDataMesh()
      self._openTimeStep(t, mesh, label, labelId)

      for name in self.vertexDataFields:
        field = self.dataProvider.getVertexField(name, fields)
        self._appendVertexField(t, field)

      for name in self.cellDataFields:
        field = self.dataProvider.getCellField(name, fields)
        self._appendCellField(t, field, label, labelId)

      self._closeTimeStep()

    self._logger.eventEnd(logEvent)
    return
      
    
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    return


  def _checkWrite(self, t):
    """
    Check if we want to write data at time t.
    """
    write = False
    
    # If first call, then _stepWrite and _tWrite are None
    if None == self._stepWrite and None == self._tWrite:
      write = True
      self._stepWrite = self._stepCur
      self._tWrite = t

    elif self.outputFreq == "skip":
      if self._stepCur > self._stepWrite + self.skip:
        write = True
        self._stepWrite = self._stepCur

    elif self.outputFreq == "time_step":
      if t >= self._tWrite + self.dt:
       write = True
       self._tWrite = t

    else:
      raise ValueError, \
            "Unknown value '%s' for output frequency." % self.outputFreq

    self._stepCur += 1
    return write


  def _verifyFields(self, available):
    """
    Verify fields for output are available.
    """
    requested = {'vertex': \
                   {'info': self.vertexInfoFields,
                    'data': self.vertexDataFields},
                 'cell': \
                   {'info': self.cellInfoFields,
                    'data': self.cellDataFields}}
    for fieldCategory in ["vertex", "cell"]:
      if not fieldCategory in available.keys():
        raise ValueError, \
            "Key '%s' not found in available fields dictionary for " \
            "object '%s'." % (fieldCategory, self.dataProvider.name)
      for dataCategory in ["info", "data"]:
        if not dataCategory in available[fieldCategory].keys():
          raise ValueError, \
              "Key '%s' not found in available fields dictionary for " \
              "object '%s'." % (fieldCategory, self.dataProvider.name)

        notavailable = []
        for name in requested[fieldCategory][dataCategory]:
          if not name in available[fieldCategory][dataCategory]:
            notavailable.append(name)
        if len(notavailable) > 0:
          msg = \
              "Requested fields not available for output.\n" \
              "Data provider: '%s'\n" \
              "Field type: '%s'\n" \
              "Data type: '%s'\n" % (self.dataProvider.name,
                                     fieldCategory, dataCategory)
          msg += "Available fields: "
          for name in available[fieldCategory][dataCategory]:
            msg += " '%s'" % name
          msg += "\n"
          msg += "Fields not available: "
          for name in notavailable:
            msg += " '%s'" % name
          raise ValueError(msg)
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Output")
    logger.initialize()

    events = ["init",
              "open",
              "close",
              "openStep",
              "closeStep",
              "writeInfo",
              "writeData"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

  def _open(self):
    raise NotImplementedError("Implement _open() in derived class.")


  def _openTimeStep(self):
    raise NotImplementedError("Implement _openTimeStep() in derived class.")


  def _appendVertexField(self):
    raise NotImplementedError("Implement _appendVertexField() in derived class.")


  def _appendCellField(self):
    raise NotImplementedError("Implement _appendCellField() in derived class.")


  def _closeTimeStep(self):
    raise NotImplementedError("Implement _closeTimeStep() in derived class.")


  def _close(self):
    raise NotImplementedError("Implement _close() in derived class.")


# End of file 
