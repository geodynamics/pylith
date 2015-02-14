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

## @file pyre/meshio/OutputManager.py
##
## @brief Python abstract base class for managing output of
## finite-element information.
##
## Factory: output_manager

from pylith.utils.PetscComponent import PetscComponent
from pylith.utils.NullComponent import NullComponent
from meshio import OutputManager as ModuleOutputManager

# OutputManager class
class OutputManager(PetscComponent, ModuleOutputManager):
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
  @li \b perf_logger Performance (memory) logger.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  from DataWriterVTK import DataWriterVTK
  writer = pyre.inventory.facility("writer", factory=DataWriterVTK, family="data_writer")
  writer.meta['tip'] = "Writer for data."

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
  
  vertexFilter = pyre.inventory.facility("vertex_filter",
                                         family="output_vertex_filter",
                                         factory=NullComponent)
  vertexFilter.meta['tip'] = "Filter for vertex data."
  
  cellFilter = pyre.inventory.facility("cell_filter",
                                       family="output_cell_filter",
                                       factory=NullComponent)
  cellFilter.meta['tip'] = "Filter for cell data."
  
  from pylith.perf.MemoryLogger import MemoryLogger
  perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                       factory=MemoryLogger)
  perfLogger.meta['tip'] = "Performance and memory logging."


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
    import weakref
    self._setupLogging()
    self.dataProvider = weakref.ref(dataProvider)
    return
  

  def verifyConfiguration(self, mesh):
    """
    Verify compatibility of configuration.
    """
    if None == self.dataProvider():
      raise ValueError("Need to set 'dataProvider' in OutputManager.")
    self._verifyFields(self.dataProvider().availableFields)

    if not "getDataMesh" in dir(self.dataProvider()):
      raise TypeError("Data provider must have a 'getDataMesh' function.")

    if len(self.vertexInfoFields) > 0 or len(self.vertexDataFields) > 0:
      if not "getVertexField" in dir(self.dataProvider()):
        raise TypeError("Data provider must have a 'getVertexField' function.")

    if len(self.cellInfoFields) > 0 or len(self.cellDataFields) > 0:
      if not "getCellField" in dir(self.dataProvider()):
        raise TypeError("Data provider must have a 'getCellField' function.")
    return


  def initialize(self, normalizer, quadrature=None):
    """
    Initialize output manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    # Nondimensionalize time step
    import weakref
    self.normalizer = normalizer
    timeScale = normalizer.timeScale()
    self.dtN = normalizer.nondimensionalize(self.dt, timeScale)

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for output is unknown."
    self.coordsys.initialize()

    if not isinstance(self.inventory.cellFilter, NullComponent):
      self.cellFilter.initialize(quadrature)
    self.writer.initialize(normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def finalize(self):
    """
    Cleanup after running problem.
    """
    if not isinstance(self.inventory.cellFilter, NullComponent):
      self.cellFilter.finalize()
    self._modelMemoryUse()
    return


  def open(self, totalTime, numTimeSteps):
    """
    Prepare for writing data.
    """
    logEvent = "%sopen" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    nsteps = self._estimateNumSteps(totalTime, numTimeSteps)

    (mesh, label, labelId) = self.dataProvider().getDataMesh()
    self._open(mesh, nsteps, label, labelId)

    self._eventLogger.eventEnd(logEvent)    
    return


  def close(self):
    """
    Perform post-write cleanup.
    """
    logEvent = "%sclose" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    self._close()

    self._eventLogger.eventEnd(logEvent)    
    return


  def writeInfo(self):
    """
    Write information fields.
    """
    logEvent = "%swriteInfo" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    if len(self.vertexInfoFields) > 0 or len(self.cellInfoFields) > 0:
      t = 0.0
      self.open(totalTime=0.0, numTimeSteps=0)
      (mesh, label, labelId) = self.dataProvider().getDataMesh()
      self._openTimeStep(t, mesh, label, labelId)

      for name in self.vertexInfoFields:
        field = self.dataProvider().getVertexField(name)
        self._appendVertexField(t, field, mesh)

      for name in self.cellInfoFields:
        field = self.dataProvider().getCellField(name)
        self._appendCellField(t, field, label, labelId)

      self._closeTimeStep()
      self._close()

    self._eventLogger.eventEnd(logEvent)
    return


  def writeData(self, t, fields):
    """
    Write fields at current time step.
    """
    logEvent = "%swriteData" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    if self._checkWrite(t) and \
           ( len(self.vertexDataFields) > 0 or \
             len(self.cellDataFields) ) > 0:

      (mesh, label, labelId) = self.dataProvider().getDataMesh()
      self._openTimeStep(t, mesh, label, labelId)

      for name in self.vertexDataFields:
        field = self.dataProvider().getVertexField(name, fields)
        self._appendVertexField(t, field, mesh)

      for name in self.cellDataFields:
        field = self.dataProvider().getCellField(name, fields)
        self._appendCellField(t, field, label, labelId)

      self._closeTimeStep()

    self._eventLogger.eventEnd(logEvent)
    return
      
    
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)

    ModuleOutputManager.coordsys(self, self.inventory.coordsys)
    ModuleOutputManager.writer(self, self.inventory.writer)
    if not isinstance(self.inventory.vertexFilter, NullComponent):
      ModuleOutputManager.vertexFilter(self, self.inventory.vertexFilter)
    if not isinstance(self.inventory.cellFilter, NullComponent):
      ModuleOutputManager.cellFilter(self, self.inventory.cellFilter)

    self.perfLogger = self.inventory.perfLogger
    return


  def _createModuleObj(self):
    """
    Create handle to C++ object.
    """
    ModuleOutputManager.__init__(self)
    return
  

  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFields('Output', self.fields())
    return


  def _estimateNumSteps(self, totalTime, numTimeSteps):
    """
    Estimate the number of time steps we expect to output.
    """
    nsteps = numTimeSteps
    if numTimeSteps > 0 and self.outputFreq == "skip" and self.skip > 0:
      nsteps = int(1 + (numTimeSteps-1) / (1+self.skip))
    elif numTimeSteps > 0 and self.outputFreq == "time_step":
      timeScale = self.normalizer.timeScale()
      totalTimeN = self.normalizer.nondimensionalize(totalTime, timeScale)
      nsteps = int(1 + totalTimeN / self.dtN)

    return nsteps


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
      if t >= self._tWrite + self.dtN:
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
            "object '%s'." % (fieldCategory, self.dataProvider().name)
      for dataCategory in ["info", "data"]:
        if not dataCategory in available[fieldCategory].keys():
          raise ValueError, \
              "Key '%s' not found in available fields dictionary for " \
              "object '%s'." % (fieldCategory, self.dataProvider().name)

        notavailable = []
        for name in requested[fieldCategory][dataCategory]:
          if not name in available[fieldCategory][dataCategory]:
            notavailable.append(name)
        if len(notavailable) > 0:
          msg = \
              "Requested fields not available for output.\n" \
              "Data provider: '%s'\n" \
              "Field type: '%s'\n" \
              "Data type: '%s'\n" % (self.dataProvider().name,
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

    self._eventLogger = logger
    return
  

  def _open(self, mesh, nsteps, label, labelId):
    """
    Call C++ open();
    """
    if label != None and labelId != None:
      ModuleOutputManager.open(self, mesh, nsteps, label, labelId)
    else:
      ModuleOutputManager.open(self, mesh, nsteps)
    return


  def _openTimeStep(self, t, mesh, label, labelId):
    """
    Call C++ openTimeStep();
    """
    if label != None and labelId != None:
      ModuleOutputManager.openTimeStep(self, t, mesh, label, labelId)
    else:
      ModuleOutputManager.openTimeStep(self, t, mesh)
    return


  def _appendVertexField(self, t, field, mesh):
    """
    Call C++ appendVertexField();
    """
    ModuleOutputManager.appendVertexField(self, t, field, mesh)
    return


  def _appendCellField(self, t, field, label, labelId):
    """
    Call C++ appendCellField();
    """
    if label != None and labelId != None:
      ModuleOutputManager.appendCellField(self, t, field, label, labelId)
    else:
      ModuleOutputManager.appendCellField(self, t, field)
    return


  def _closeTimeStep(self):
    """
    Call C++ closeTimeStep().
    """
    ModuleOutputManager.closeTimeStep(self)
    return


  def _close(self):
    """
    Call C++ close().
    """
    ModuleOutputManager.close(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputManager()


# End of file 
