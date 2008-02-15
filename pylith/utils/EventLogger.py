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

## @file pylith/utils/EventLogger.py
##
## @brief Python object for managing event logging using PETSc.
##
## Each logger object manages the events for a single "logging class".

# EventLogger class
class EventLogger(object):
  """
  Python object for managing event logging using PETSc.

  Each logger object manages the events for a single 'logging class'.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    self.cppHandle = None
    self._createCppHandle()
    self.events = {} # dict of events with counts for current logging.
    return


  def setClassName(self, name):
    """
    Set name of logging class.
    """
    self._createCppHandle()
    self.cppHandle.className = name
    return


  def getClassName(self):
    """
    Set name of logging class.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.className


  def initialize(self):
    """
    Setup logging class.
    """    
    assert(None != self.cppHandle)
    self.cppHandle.initialize()
    return


  def registerEvent(self, name):
    """
    Register event.
    """
    self.events[name] = 0 # Set log count to 0
    assert(None != self.cppHandle)
    return self.cppHandle.registerEvent(name)


  def eventId(self, name):
    """
    Get event identifier.
    """    
    assert(None != self.cppHandle)
    return self.cppHandle.eventId(name)


  def eventBegin(self, name):
    """
    Log event begin.
    """
    if self.events[name] == 0: # prevent double counting
      assert(None != self.cppHandle)
      self.cppHandle.eventBegin(self.cppHandle.eventId(name))
    self.events[name] += 1
    return


  def eventEnd(self, name):
    """
    Log event end.
    """
    if self.events[name] > 0:
      self.events[name] -= 1
    if 0 == self.events[name]: # prevent double counting
      assert(None != self.cppHandle)
      self.cppHandle.eventEnd(self.cppHandle.eventId(name))
    
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.utils.utils as bindings
      self.cppHandle = bindings.EventLogger()
  
  
# End of file 
