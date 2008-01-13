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

## @file pyre/meshio/SolutionIO.py
##
## @brief Python abstract base class for I/O of the finite-element
## solution.
##
## Factory: solution_io

from pyre.components.Component import Component

# SolutionIO class
class SolutionIO(Component):
  """
  Python abstract base class for finite-element mesh I/O.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing SolutionIOVTK facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SolutionIOVTK facilities and properties.
    ##
    ## \b Properties
    ## @li \b output_freq Flag indicating whether to use 'time_step' or 'skip'
    ##   to set frequency of solution output.
    ## @li \b time_step Time step between solution output.
    ## @li \b skip Number of time steps to skip between solution output.
    ##
    ## \b Facilities
    ## @li \b coordsys Coordinate system for output.

    import pyre.inventory

    outputFreq = pyre.inventory.str("output_freq", default="skip",
             validator=pyre.inventory.choice(["skip", "time_step"]))
    outputFreq.meta['tip'] = "Flag indicating whether to use 'time_step' " \
                             "or 'skip' to set frequency of solution output."

    from pyre.units.time import s
    dt = pyre.inventory.dimensional("time_step", default=1.0*s)
    dt.meta['tip'] = "Time step between solution output."

    skip = pyre.inventory.int("skip", default=0,
                              validator=pyre.inventory.greaterEqual(0))
    skip.meta['tip'] = "Number of time steps to skip between solution output."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system for output."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solutionio"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution_io")
    self.cppHandle = None
    self.coordsys = None
    self.mesh = None
    self.t = None
    self.istep = None
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def open(self, mesh):
    """
    Open files for solution.
    """
    self._setupLogging()
    logEvent = "%sopen" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    
    
    self._info.log("Opening files for output of solution.")
    self.mesh = mesh

    # Set flags
    self._sync()

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for output is unknown."
    self.coordsys.initialize()

    assert(self.cppHandle != None)
    self.cppHandle.open(mesh.cppHandle)

    self._logger.eventEnd(logEvent)    
    return


  def close(self):
    """
    Close files for solution.
    """
    logEvent = "%sclose" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    
    self._info.log("Closing files for output of solution.")

    # Set flags
    self._sync()

    assert(self.cppHandle != None)
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
    if self.istep == None or not "value" in dir(self.t):
      write = True
    elif self.outputFreq == "skip":
      if istep > self.istep + self.skip:
        write = True
    elif t >= self.t + self.dt:
      write = True
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


  def writeVertexField(self, t, istep, name, field):
    """
    Write field over vertices at time t to file.
    """
    logEvent = "%swriteVertex" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    if self.writeFlag:
      self._info.log("Writing solution field '%s'." % name)
      assert(self.cppHandle != None)
      assert(self.mesh.cppHandle != None)
      self.cppHandle.writeVertexField(t.value, name, field,
                                      self.mesh.cppHandle)
      self.istep = istep
      self.t = t

    self._logger.eventEnd(logEvent)
    return


  def writeCellField(self, t, istep, name, field):
    """
    Write field over cells at time t to file.
    """
    logEvent = "%swriteCell" % self._loggingPrefix
    self._logger.eventBegin(logEvent)    

    if self.writeFlag:
      self._info.log("Writing solution field '%s'." % name)
      assert(self.cppHandle != None)
      assert(self.mesh.cppHandle != None)
      self.cppHandle.writeCellField(t.value, name, field, 
                                    self.mesh.cppHandle)
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
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    self._createCppHandle()
    self.cppHandle.coordsys = self.coordsys.cppHandle
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createCppHandle() in " \
                              "derived class.")
  
  
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
              "openStep",
              "closeStep",
              "writeVertex",
              "writeCell"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# End of file
