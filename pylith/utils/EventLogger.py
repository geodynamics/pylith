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
    self.event = None
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
    if self.event != name: # prevent double logging
      assert(None != self.cppHandle)
      self.cppHandle.eventBegin(self.cppHandle.eventId(name))
      self.event = name
    return


  def eventEnd(self, name):
    """
    Log event end.
    """
    if None != self.event: # prevent double logging
      assert(None != self.cppHandle)
      self.cppHandle.eventEnd(self.cppHandle.eventId(name))
      self.event = None
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
